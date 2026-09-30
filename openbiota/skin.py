"""Gut-microbiome research patterns for inflammatory skin disease — spec v4.2.

Nine numeric panels across six clinical families, three registered tasks that
deliberately produce no number, and the counterevidence that stops any of them
being read as a diagnosis.

What this module is not
-----------------------
There is no established, uniquely diagnostic stool signature for any of these
conditions. Every panel here is a *source-cohort association*: a list of
organisms that differed between one study's cases and one study's controls.
A high index means this sample's community moves in the direction that study
reported. It does not mean the person has the disease, and a low index does
not mean they do not. Healthy-versus-case discrimination in a single centre
is not specificity; that would need testing against other skin diseases, IBD,
medications and populations, which none of these sources did.

Several things follow from that, and are enforced rather than documented:

* Panels in one family are never averaged into a family score. They are shown
  side by side with their own contrasts and evidence tiers, because a pooled
  psoriasis/PsA contrast cannot supply a PsA headline (§5.1, E04).
* A feature appearing in several panels is one observation, not several. The
  taxonomy and function panels from the same participants are one evidence
  group (:attr:`Panel.same_participants_as`).
* Published AUCs belong to the published classifiers. None of them is
  attached to the index computed here, which is a different estimator on
  different features (E01, E04, E08).
* SGB4348 is not psoriasis-specific: the same marker is raised in pemphigus
  foliaceus (E03). The counterassociation travels with the panel.
* Seborrheic dermatitis gets a page and no score. Culture counts and
  Mendelian-randomisation directions are real evidence about different
  quantities, and neither converts into a DNA relative-abundance coefficient.

The arithmetic, missingness policy and reference fitting are
:mod:`openbiota.dirindex`, shared with the chronic-hives panels.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any, Final

from openbiota.dirindex import (
    DirIndexError,
    FrozenReference,
    Panel,
    PanelFeature,
    PanelResult,
    freeze_reference,
    score_panel,
)

RELEASE: Final = "4.2"

#: Score semantics per measurement mode.
NATIVE: Final = "unvalidated_directional_pattern_index"
TRANSPORTED: Final = "unvalidated_16S_to_shotgun_directional_pattern_index"
KO: Final = "unvalidated_directional_ko_pattern_index"

#: Engineering QC floor for an ordinary panel (§4.3). Pairs, not reads.
MIN_USABLE_READ_PAIRS: Final = 500_000

#: Recent-antibiotic exclusion window for participant-facing numbers (§4.3).
ANTIBIOTIC_EXCLUSION_DAYS: Final = 90


# --------------------------------------------------------------------------- #
# families
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class Family:
    """A clinical phenotype, which may carry several study panels."""

    family_id: str
    label: str
    #: Always names the condition in words a reader recognises.
    plain_label: str
    summary: str
    #: What a high index does and does not mean here.
    meaning: str

    def to_json(self) -> dict[str, Any]:
        return {
            "family_id": self.family_id, "label": self.label,
            "plain_label": self.plain_label, "summary": self.summary,
            "result_meaning": self.meaning,
        }


FAMILIES: Final[dict[str, Family]] = {
    "plaque_psoriasis": Family(
        family_id="plaque_psoriasis",
        label="Plaque psoriasis",
        plain_label="plaque psoriasis (raised scaly skin plaques)",
        summary=(
            "Plaque psoriasis is a chronic immune-mediated skin disease with raised, scaly plaques, "
            "commonly on elbows, knees, scalp and trunk. Four separate stool-shotgun studies are carried "
            "here as separate readings rather than one score, because they disagree: Deng 2026 "
            "nominates a single SGB, Chang 2022 and Sá 2026 report largely non-overlapping species, "
            "and Bacteroides coprocola is enriched in cases in one study and in controls in another."
        ),
        meaning=(
            "Research pattern concordance. A high reading does not establish psoriasis, and the marker "
            "with the strongest single-feature result here is also raised in a different blistering "
            "skin disease."
        ),
    ),
    "psoriatic_disease_shared": Family(
        family_id="psoriatic_disease_shared",
        label="Psoriatic disease (skin and joint forms pooled)",
        plain_label="psoriatic disease, skin and joint forms pooled",
        summary=(
            "Xiao 2024 compared 44 people with skin-only plaque psoriasis and 26 with psoriatic "
            "arthritis, pooled, against 25 controls. The three species below come from that pooled "
            "comparison. No species difference separated psoriatic arthritis from skin-only psoriasis "
            "in that cohort, so this reading cannot tell the two apart."
        ),
        meaning=(
            "Resemblance to a pooled psoriatic-disease contrast. It is shown under both the psoriasis "
            "and the psoriatic-arthritis pages, and it is not a subtype test on either."
        ),
    ),
    "psoriatic_arthritis": Family(
        family_id="psoriatic_arthritis",
        label="Psoriatic arthritis",
        plain_label="psoriatic arthritis (psoriasis with joint inflammation)",
        summary=(
            "Psoriatic arthritis is inflammatory arthritis occurring with psoriasis. The panel here is "
            "nine microbial gene families (KEGG orthologues) reported by Liu 2024 in 20 treatment-naive "
            "patients against 10 controls. Ten controls do not meet this build's reference minimum, so "
            "the panel scores only against an independently compatible functional reference, or not at all."
        ),
        meaning=(
            "Research resemblance on directly measured gene-family abundances. It is not a validated "
            "test distinguishing psoriatic arthritis from skin-only psoriasis, and it says nothing "
            "about joints."
        ),
    ),
    "adult_atopic_dermatitis": Family(
        family_id="adult_atopic_dermatitis",
        label="Atopic dermatitis, adult",
        plain_label="atopic dermatitis in adults (eczema)",
        summary=(
            "Atopic dermatitis, commonly called eczema, is a chronic itchy inflammatory skin disease. "
            "Wang 2023 studied 104 adult cases and 130 controls by 16S amplicon sequencing, which names "
            "organisms differently from the shotgun sequencing used here. Thirteen slots are declared; "
            "five are source groups from a 16S reference database with no exact shotgun equivalent, and "
            "those stay in the coverage denominator as unavailable rather than being approximated."
        ),
        meaning=(
            "Exploratory translated pattern. It inherits none of the source study's sensitivity or "
            "specificity, because it is not the same measurement."
        ),
    ),
    "nonsegmental_vitiligo": Family(
        family_id="nonsegmental_vitiligo",
        label="Non-segmental vitiligo",
        plain_label="non-segmental vitiligo (patches of skin losing pigment)",
        summary=(
            "Vitiligo is loss of skin pigment in patches; the non-segmental form is the commoner, "
            "usually symmetrical type. Two panels are carried separately. Luan 2023 is a single matched "
            "cohort of 25 cases and 25 controls with native shotgun sequencing. Ju 2025 compared 10 "
            "actively spreading cases against 20 controls selected out of a different study, so its "
            "cases and controls were never processed together."
        ),
        meaning=(
            "Cohort-specific resemblance. No diagnostic probability, and the two panels are not "
            "replications of each other — they share almost no features."
        ),
    ),
    "hidradenitis_suppurativa": Family(
        family_id="hidradenitis_suppurativa",
        label="Hidradenitis suppurativa",
        plain_label="hidradenitis suppurativa (recurrent painful skin abscesses in folds)",
        summary=(
            "Hidradenitis suppurativa causes recurrent painful nodules and abscesses in skin folds. The "
            "stool evidence is thin: one 16S study of 15 cases and 15 controls supplies a single genus "
            "at p=0.046, and a separate shotgun study found no significant difference at all between "
            "hidradenitis cases and controls. Both are shown."
        ),
        meaning=(
            "A one-feature translated reading. With a single slot there is no pattern to speak of, and "
            "the negative shotgun comparison is the more informative result on this page."
        ),
    ),
    "seborrheic_dermatitis": Family(
        family_id="seborrheic_dermatitis",
        label="Seborrheic dermatitis",
        plain_label="seborrheic dermatitis (scaly, greasy scalp and facial rash; dandruff at its mildest)",
        summary=(
            "Seborrheic dermatitis is a scaly, greasy rash of the scalp, face and upper trunk; dandruff "
            "is its mildest form. Gut involvement has been reported, but the available evidence is stool "
            "bacterial culture with no per-subject values, genetically instrumented associations rather "
            "than measured cases, a clinical series with no sequencing, and scalp-skin studies. None of "
            "those defines a stool-DNA signature."
        ),
        meaning=(
            "No score is computed, and that is not a negative test. It means no adequately specified "
            "stool-DNA signature was identified for this implementation."
        ),
    ),
    "infant_atopic_dermatitis": Family(
        family_id="infant_atopic_dermatitis",
        label="Atopic dermatitis, infant (strain-level task)",
        plain_label="atopic dermatitis in infants (infant eczema), strain-level",
        summary=(
            "Seong 2025 found no significant species-abundance differences at six months. The signal was "
            "a subclade distinction inside Bifidobacterium longum subspecies infantis, which needs "
            "strain-level reconstruction against released assemblies, not a species abundance."
        ),
        meaning=(
            "Registered and pending. It is a six-month-infant finding and is never applied to an adult "
            "or to a child aged 8 to 18."
        ),
    ),
}


# --------------------------------------------------------------------------- #
# panels (§3.1 canonical seed)
# --------------------------------------------------------------------------- #

_SGB4348_REASON: Final = (
    "Requires the same biological SGB definition as its sources. A database without that SGB "
    "returns unavailable: the same number in another namespace is a different organism, and no "
    "phylum, family or unclassified bucket substitutes for it."
)

PSO_DENG2026: Final = Panel(
    panel_id="PSO_DENG2026_SGB4348_V1",
    revision="4.2.0",
    label="Plaque psoriasis — Deng 2026 single-marker SGB4348",
    family="plaque_psoriasis",
    contrast="psoriasis_vs_controls",
    source_id="DENG_2026",
    source_ids=("E01", "E02"),
    counterevidence_ids=("E03", "E05"),
    source_assay="shotgun_metagenomics",
    scoring_assay="shotgun_metagenomics",
    measurement="shotgun_taxonomic_fraction",
    score_semantics=NATIVE,
    evidence_tier="same_center_tested_candidate",
    discovery_sample_count=143,
    percentile_note=(
        "The source reported a single-feature ROC AUC of 0.84 (95% CI 0.74-0.94) on a same-centre "
        "test group. That figure belongs to the source's own threshold on its own cohort. It is not "
        "the accuracy of this index and is never attached to an individual result."
    ),
    features=(
        PanelFeature(
            "s__GGB51647_SGB4348", 1, rank="sgb",
            measured_as="GGB51647_SGB4348",
            mapping_reason=_SGB4348_REASON,
        ),
    ),
)

PSO_CHANG2022: Final = Panel(
    panel_id="PSO_CHANG2022_TAX_V1",
    revision="4.2.0",
    label="Plaque psoriasis — Chang 2022 seven-taxon panel",
    family="plaque_psoriasis",
    contrast="psoriasis_vs_healthy",
    source_id="CHANG_2022",
    source_ids=("E14",),
    counterevidence_ids=("E04",),
    source_assay="shotgun_metagenomics",
    scoring_assay="shotgun_metagenomics",
    measurement="shotgun_taxonomic_fraction",
    score_semantics=NATIVE,
    evidence_tier="single_cohort_selected_subset",
    discovery_sample_count=48,
    features=(
        PanelFeature(
            "Megasphaera_unclassified", 1, rank="source_group", q_value=0.004186369,
            mapping_reason=(
                "The source's unclassified Megasphaera bucket, not the genus total and not "
                "M. elsdenii. Without that exact bucket in the production catalogue the slot is "
                "unavailable."
            ),
        ),
        PanelFeature("Dialister_invisus", 1, q_value=0.001738535),
        PanelFeature(
            "Bacteroides_vulgatus", 1, q_value=0.01176506,
            synonyms=("Phocaeicola_vulgatus",),
            mapping_reason=(
                "Bacteroides vulgatus was reclassified to Phocaeicola vulgatus (García-López 2019): "
                "one organism, two names. The pinned catalogue is checked for both names and the one "
                "it carries is used; if it carried both, the slot would be withheld rather than summed."
            ),
        ),
        PanelFeature("Phascolarctobacterium_succinatutens", -1, q_value=0.040279303),
        PanelFeature("Bifidobacterium_pseudocatenulatum", -1, q_value=0.002393849),
        PanelFeature(
            "Coprococcus_sp_ART55_1", -1, rank="source_species", q_value=0.000961438,
            mapping_reason=(
                "An exact source concept, not the genus Coprococcus and not any named Coprococcus "
                "species. Unavailable unless the catalogue carries this specific entry."
            ),
        ),
        # Direction conflict with E04, which found this organism enriched in
        # controls. Kept as the source reported it, and flagged on the page.
        PanelFeature("Bacteroides_coprocola", 1, q_value=1.33958e-10),
    ),
)

PSO_SA2026: Final = Panel(
    panel_id="PSO_SA2026_TAX_V1",
    revision="4.2.0",
    label="Severe plaque psoriasis — Sá 2026 three-taxon panel",
    family="plaque_psoriasis",
    contrast="severe_psoriasis_vs_healthy",
    source_id="SA_2026",
    source_ids=("E03",),
    source_assay="shotgun_metagenomics",
    scoring_assay="shotgun_metagenomics",
    measurement="shotgun_taxonomic_fraction",
    score_semantics=NATIVE,
    evidence_tier="small_confounded_cohort",
    discovery_sample_count=34,
    percentile_note=(
        "Twenty-four psoriasis participants, of whom 18 were receiving immunosuppression and 5 had "
        "parasitic infection, against 10 controls. Treatment and infection are not separable from "
        "disease in this contrast."
    ),
    features=(
        PanelFeature(
            "Enterocloster", -1, rank="genus",
            aggregate_of=(
                "Clostridium_bolteae|Enterocloster_bolteae",
                "Clostridium_bolteae_CAG_59",
                "Clostridium_clostridioforme|Enterocloster_clostridioformis",
                "Clostridium_aldenense|Enterocloster_aldenensis",
                "Clostridium_citroniae|Enterocloster_citroniae",
                "Clostridium_lavalense|Enterocloster_lavalensis",
                "Clostridium_asparagiforme|Enterocloster_asparagiformis",
            ),
            mapping_reason=(
                "Enterocloster was carved out of Clostridium in 2020 (Haas and Blanchard) from six "
                "named species. The genus total is the sum of exactly those species under whichever "
                "name the pinned catalogue carries; a catalogue carrying both names of one species "
                "withholds the slot rather than double counting. No other Clostridium is included."
            ),
        ),
        PanelFeature("Phascolarctobacterium_faecium", -1),
        PanelFeature(
            "Bacteroides_finegoldii", -1,
            mapping_reason=(
                "Bacteroides finegoldii, which is a different organism from Alistipes finegoldii in "
                "the pooled psoriatic-disease panel. Same epithet, different genus; no fuzzy match."
            ),
        ),
    ),
)

PSORIATIC_SHARED_XIAO2024: Final = Panel(
    panel_id="PSORIATIC_SHARED_XIAO2024_TAX_V1",
    revision="4.2.0",
    label="Psoriatic disease, skin and joint pooled — Xiao 2024 three species",
    family="psoriatic_disease_shared",
    contrast="pooled_plaque_psoriasis_and_psoriatic_arthritis_vs_healthy",
    source_id="XIAO_2024",
    source_ids=("E04",),
    source_assay="shotgun_metagenomics",
    scoring_assay="shotgun_metagenomics",
    measurement="shotgun_taxonomic_fraction",
    score_semantics=NATIVE,
    evidence_tier="shared_psoriasis_signal",
    discovery_sample_count=95,
    min_age_years=18.0,
    max_age_years=65.0,
    percentile_note=(
        "The q-values here belong to the pooled psoriatic-disease-versus-healthy contrast and to "
        "nothing else. The source reported roughly 0.74-0.77 AUC for Eubacterium rectale abundance "
        "alone; that is one feature's figure, not this three-species index's, and subgroup-specific "
        "estimates for the two Alistipes species are unknown."
    ),
    features=(
        PanelFeature(
            "Eubacterium_rectale", -1, q_value=0.01,
            synonyms=("Agathobacter_rectalis",),
            mapping_reason=(
                "Eubacterium rectale was transferred to Agathobacter rectalis (Rosero 2016): one "
                "organism, two names. The pinned catalogue is checked for both and the one it carries "
                "is used; a catalogue carrying both would withhold the slot rather than sum them."
            ),
        ),
        PanelFeature("Alistipes_finegoldii", -1, q_value=0.001),
        PanelFeature("Alistipes_shahii", -1, q_value=0.001),
    ),
)

PSA_LIU2024_KO: Final = Panel(
    panel_id="PSA_LIU2024_KO_V1",
    revision="4.2.0",
    label="Psoriatic arthritis — Liu 2024 nine gene families",
    family="psoriatic_arthritis",
    contrast="treatment_naive_psoriatic_arthritis_vs_healthy",
    source_id="LIU_2024",
    source_ids=("E15",),
    source_assay="shotgun_metagenomics",
    scoring_assay="shotgun_ko_annotation",
    measurement="shotgun_ko_fraction",
    score_semantics=KO,
    evidence_tier="single_cohort_nominal",
    discovery_sample_count=30,
    min_age_years=18.0,
    max_age_years=60.0,
    percentile_note=(
        "Selected by LEfSe at p<0.05 with LDA>2 in 20 patients and 10 controls, with no external "
        "validation and no recovered public accession. Ten discovery controls do not meet this "
        "build's reference minimum of 20 per feature."
    ),
    features=(
        PanelFeature("K02004", 1, rank="ko"),
        PanelFeature("K01190", 1, rank="ko"),
        PanelFeature("K05349", 1, rank="ko"),
        PanelFeature("K01897", 1, rank="ko"),
        PanelFeature("K01187", 1, rank="ko"),
        PanelFeature(
            "K07114", 1, rank="ko",
            mapping_reason=(
                "K07114 as printed in the source's results. A conflicting spelling K07714 appears in "
                "its discussion and is not imported. This is a microbial calcium-activated-chloride-"
                "channel homologue annotation: it is not a human chloride channel."
            ),
        ),
        PanelFeature("K06131", -1, rank="ko"),
        PanelFeature("K07729", -1, rank="ko"),
        PanelFeature("K01005", -1, rank="ko"),
    ),
)

ATD_ADULT_WANG2023: Final = Panel(
    panel_id="ATD_ADULT_WANG2023_TRANSLATED_V1",
    revision="4.2.0",
    label="Adult atopic dermatitis (eczema) — Wang 2023 thirteen slots, translated from 16S",
    family="adult_atopic_dermatitis",
    contrast="adult_atopic_dermatitis_vs_healthy",
    source_id="WANG_2023",
    source_ids=("E06",),
    counterevidence_ids=("E07",),
    source_assay="16S_amplicon_V4_V5",
    scoring_assay="shotgun_metagenomics",
    measurement="shotgun_genus_fraction_translated_from_16s",
    score_semantics=TRANSPORTED,
    assay_transport="unvalidated_16s_to_shotgun",
    evidence_tier="assay_transported_exploratory",
    discovery_sample_count=234,
    min_age_years=18.0,
    max_age_years=65.0,
    percentile_note=(
        "Selected by LEfSe at p<0.05 with LDA>2. Only Clostridium sensu stricto 1 was also "
        "significant under ANCOM (W=218). Neither Shannon diversity nor the Firmicutes-to-"
        "Bacteroidetes ratio differed, so neither appears here."
    ),
    features=(
        PanelFeature("Blautia", 1, rank="genus"),
        PanelFeature("Butyricicoccus", 1, rank="genus"),
        PanelFeature("Lachnoclostridium", 1, rank="genus"),
        PanelFeature(
            "Eubacterium_hallii_group", 1, rank="source_group",
            mapping_reason=(
                "A SILVA138 group, not a genus and not a species. Its membership must be resolved to "
                "a frozen non-overlapping set of catalogue entries before it can be summed; until "
                "then the slot is unavailable and stays in the denominator."
            ),
        ),
        PanelFeature("Erysipelatoclostridium", 1, rank="genus"),
        PanelFeature("Megasphaera", 1, rank="genus"),
        PanelFeature("Oscillibacter", 1, rank="genus"),
        PanelFeature("Flavonifractor", 1, rank="genus"),
        PanelFeature(
            "unclassified_Oscillospiraceae", 1, rank="source_group",
            mapping_reason=(
                "The source's unclassified bucket within this family. It is not the family total: "
                "summing the whole family would count organisms the source resolved separately."
            ),
        ),
        PanelFeature("Romboutsia", -1, rank="genus"),
        PanelFeature(
            "Clostridium_sensu_stricto_1", -1, rank="source_group",
            mapping_reason=(
                "A SILVA138 group, and the one finding here that also survived ANCOM. It is never "
                "substituted with the genus Clostridium, which is far broader."
            ),
        ),
        PanelFeature(
            "unclassified_Butyricicoccaceae", -1, rank="source_group",
            mapping_reason="The source's unclassified bucket within this family, not the family total.",
        ),
        PanelFeature(
            "unclassified_Erysipelotrichaceae", -1, rank="source_group",
            mapping_reason="The source's unclassified bucket within this family, not the family total.",
        ),
    ),
)

VIT_NONSEG_LUAN2023: Final = Panel(
    panel_id="VIT_NONSEG_LUAN2023_TAX_V1",
    revision="4.2.0",
    label="Advanced non-segmental vitiligo — Luan 2023 eleven species",
    family="nonsegmental_vitiligo",
    contrast="advanced_nonsegmental_vitiligo_vs_healthy",
    source_id="LUAN_2023",
    source_ids=("E08",),
    source_assay="shotgun_metagenomics",
    scoring_assay="shotgun_metagenomics",
    measurement="shotgun_taxonomic_fraction",
    score_semantics=NATIVE,
    evidence_tier="single_cohort_selected_subset",
    discovery_sample_count=50,
    percentile_note=(
        "The source's published AUC of 0.786 belongs to a 28-feature random forest under nested "
        "cross-validation, which is a different estimator on a different feature set from the "
        "eleven-slot index computed here. Beta diversity did not differ (p=0.495). Individual "
        "q-values were not published per feature and are not invented."
    ),
    features=(
        PanelFeature(
            "Bacteroides_fragilis", 1,
            mapping_reason=(
                "Species-level relative abundance. It is not an enterotoxigenic strain call: the "
                "toxin gene is not measured by a species abundance."
            ),
        ),
        PanelFeature("Massilioclostridium_coli", 1),
        PanelFeature(
            "Lachnospiraceae_bacterium_BX3", 1, rank="source_species",
            mapping_reason="An exact source concept, not the family Lachnospiraceae.",
        ),
        PanelFeature(
            "TM7_phylum_sp_oral_taxon_348", 1, rank="source_species",
            mapping_reason="An exact source concept, not the TM7 phylum as a whole.",
        ),
        PanelFeature(
            "Prevotella_copri_clade_C", -1, rank="source_clade",
            mapping_reason=(
                "One clade of Prevotella copri, not all P. copri and not an unverified Segatella "
                "relabelling. Without clade resolution in the catalogue the slot is unavailable."
            ),
        ),
        PanelFeature("Dorea_longicatena", -1),
        PanelFeature("Allisonella_histaminiformans", -1),
        PanelFeature("Coprobacter_secundus", -1),
        PanelFeature("Coprococcus_comes", -1),
        PanelFeature("Sellimonas_intestinalis", -1),
        PanelFeature("Bacteroides_bouchesdurhonensis", -1),
    ),
)

VIT_ACTIVE_JU2025: Final = Panel(
    panel_id="VIT_ACTIVE_JU2025_TAX_V1",
    revision="4.2.0",
    label="Actively spreading vitiligo — Ju 2025 four species, exploratory",
    family="nonsegmental_vitiligo",
    contrast="active_spreading_vitiligo_vs_external_controls",
    source_id="JU_2025",
    source_ids=("E09",),
    source_assay="shotgun_metagenomics",
    scoring_assay="shotgun_metagenomics",
    measurement="shotgun_taxonomic_fraction",
    score_semantics=NATIVE,
    evidence_tier="cross_study_control_confounding",
    discovery_sample_count=30,
    percentile_note=(
        "Ten cases from one project and twenty controls selected out of a different project of 37 "
        "people. Cases and controls were never processed together, so study batch and disease are "
        "not separable. Alpha diversity did not differ significantly."
    ),
    features=(
        PanelFeature(
            "Faecalibacterium_prausnitzii", -1,
            mapping_reason=(
                "Modern taxonomy splits historical F. prausnitzii, and F. duncaniae is one of the "
                "organisms split out of it. Both slots are scored only when the catalogue defines "
                "them disjointly; otherwise the pair is quarantined rather than summed or "
                "substituted, because an umbrella species and one of its children are not two "
                "independent observations."
            ),
        ),
        PanelFeature(
            "Faecalibacterium_duncaniae", -1,
            mapping_reason=(
                "Split out of historical F. prausnitzii. Quarantined together with it unless the "
                "catalogue proves the two inputs disjoint."
            ),
        ),
        PanelFeature("Megamonas_funiformis", -1),
        PanelFeature("Bifidobacterium_bifidum", 1),
    ),
)

HS_OGUT2022: Final = Panel(
    panel_id="HS_OGUT2022_TRANSLATED_V1",
    revision="4.2.0",
    label="Hidradenitis suppurativa — Ogut 2022 single genus, translated from 16S",
    family="hidradenitis_suppurativa",
    contrast="hidradenitis_suppurativa_vs_healthy",
    source_id="OGUT_2022",
    source_ids=("E10",),
    counterevidence_ids=("E03",),
    source_assay="16S_amplicon_V3_V4",
    scoring_assay="shotgun_metagenomics",
    measurement="shotgun_genus_fraction_translated_from_16s",
    score_semantics=TRANSPORTED,
    assay_transport="unvalidated_16s_to_shotgun",
    evidence_tier="single_feature_nominal",
    discovery_sample_count=30,
    percentile_note=(
        "One genus at p=0.046 in 15 cases and 15 controls, by 16S amplicon, with no verified raw "
        "accession. A separate shotgun study (Sá 2026) found no significant taxonomic difference "
        "between hidradenitis cases and controls at all."
    ),
    features=(
        PanelFeature("Fusicatenibacter", -1, rank="genus", q_value=None),
    ),
)

#: KEGG:00072 decreased in psoriatic arthritis versus skin-only psoriasis
#: (OR 0.22, q=0.044). A subtype contrast, not a healthy-screening signal, so
#: it is registered as its own panel and never folded into either family score.
PSA_VS_PSO_KEGG00072: Final = Panel(
    panel_id="PSA_VS_PSO_XIAO2024_KEGG00072_V1",
    revision="4.2.0",
    label="Psoriatic arthritis versus skin-only psoriasis — KEGG 00072 capacity",
    family="psoriatic_arthritis",
    contrast="psoriatic_arthritis_vs_skin_only_psoriasis",
    source_id="XIAO_2024",
    source_ids=("E04",),
    source_assay="shotgun_metagenomics",
    scoring_assay="shotgun_pathway_annotation",
    measurement="shotgun_ko_fraction",
    score_semantics=KO,
    evidence_tier="subtype_contrast_capacity_observation",
    discovery_sample_count=70,
    min_age_years=18.0,
    max_age_years=65.0,
    same_participants_as=("PSORIATIC_SHARED_XIAO2024_TAX_V1",),
    percentile_note=(
        "This contrast compares two groups who both have psoriatic disease. It answers nothing about "
        "a healthy person and is never inferred from the three taxa in the pooled panel from the same "
        "participants. Butanoate metabolism in the same source was not significant (q=0.157)."
    ),
    features=(
        PanelFeature("KEGG:00072", -1, rank="pathway"),
    ),
)

#: Every numeric panel in this release, in report order.
PANELS: Final[tuple[Panel, ...]] = (
    PSO_DENG2026,
    PSO_CHANG2022,
    PSO_SA2026,
    PSORIATIC_SHARED_XIAO2024,
    PSA_LIU2024_KO,
    PSA_VS_PSO_KEGG00072,
    ATD_ADULT_WANG2023,
    VIT_NONSEG_LUAN2023,
    VIT_ACTIVE_JU2025,
    HS_OGUT2022,
)

#: The nine panels the specification counts as the numeric seed inventory
#: (§3.1). The subtype-contrast panel above is an additional registered
#: observation, not part of that count.
SEED_PANELS: Final[tuple[Panel, ...]] = tuple(p for p in PANELS if p is not PSA_VS_PSO_KEGG00072)
SEED_SLOT_COUNT: Final = 52

BY_ID: Final[dict[str, Panel]] = {p.panel_id: p for p in PANELS}


# --------------------------------------------------------------------------- #
# tasks that deliberately produce no number (§3.1 unavailable_tasks)
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class UnavailableTask:
    """A registered task with no defensible number at this cutoff.

    Visible, with a truthful reason. An unavailable score is not a negative
    result and not clearance: it means the measurement or the reference does
    not exist yet.
    """

    task_id: str
    family: str
    state: str
    source_ids: tuple[str, ...]
    #: Exact participant-facing text where the specification fixes it.
    statement: str
    #: What would have to exist for this task to produce a number.
    promotion_requirement: str
    #: Plain-language heading for the report.
    label: str = ""

    def to_json(self) -> dict[str, Any]:
        return {
            "panel_id": self.task_id,
            "label": self.label or self.task_id.replace("_", " "),
            "profile_family": self.family,
            "status": self.state,
            "index": None,
            "index_label": "Research directional pattern index",
            "usable_count": 0,
            "panel_count": 0,
            "coverage": 0.0,
            "interval_coverage": 0.0,
            "feature_mask": [],
            "full_panel_bounds": [0.0, 100.0],
            "bounds_type": "coverage_and_censoring",
            "control_percentile": None,
            "case_percentile": None,
            "disease_probability": None,
            "clinical_classification": None,
            "disease_specificity": "not_established",
            "source_ids": list(self.source_ids),
            "reason_codes": [self.state],
            "statement": self.statement,
            "promotion_requirement": self.promotion_requirement,
            "feature_results": [],
        }


SEBD_STATEMENT: Final = (
    "Gut involvement has been reported, but no adequately specified stool-DNA signature was "
    "identified for this implementation. This result cannot assess whether you have seborrheic "
    "dermatitis."
)

UNAVAILABLE_TASKS: Final[tuple[UnavailableTask, ...]] = (
    UnavailableTask(
        task_id="SEBD_STOOL_RESEARCH_V1",
        label="Seborrheic dermatitis — the evidence, and why there is no score",
        family="seborrheic_dermatitis",
        state="insufficient_stool_signature",
        source_ids=("E11", "E12", "E13"),
        statement=SEBD_STATEMENT,
        promotion_requirement=(
            "A measured seborrheic-dermatitis case-control stool dataset with exact phenotype, body "
            "site, feature definitions, assay, clinical metadata and a reproducible analysis, tested "
            "against psoriasis and sebopsoriasis, atopic dermatitis, rosacea and healthy controls. A "
            "culture-based module would additionally need quantitative culture inputs and reference "
            "distributions; it is not built by converting published percentages into colony counts."
        ),
    ),
    UnavailableTask(
        task_id="ATD_INFANT_SEONG2025_STRAIN_V1",
        label="Infant atopic dermatitis (eczema) — Seong 2025 strain-level finding, not scored",
        family="infant_atopic_dermatitis",
        state="strain_bundle_pending",
        source_ids=("E07",),
        statement=(
            "The infant finding here is a distinction between two subclades inside one bacterial "
            "subspecies, which this pipeline cannot yet resolve. No species or genus measurement "
            "substitutes for it, and the finding does not transfer to an adult or to a child aged "
            "8 to 18."
        ),
        promotion_requirement=(
            "The exact assemblies and marker reference, source subclade labels, sequence quality and "
            "unique-marker alignment criteria, strain coverage and mixture behaviour, and independent "
            "age-matched validation. No read-depth threshold is asserted for strain reconstruction "
            "until such a validation bundle exists."
        ),
    ),
    UnavailableTask(
        task_id="PSO_DENG2026_ORIGINAL_MODEL",
        label="Plaque psoriasis — Deng 2026 original classifier, not reproduced",
        family="plaque_psoriasis",
        state="blocked_unresolved_preprocessing_and_artifacts",
        source_ids=("E01", "E02"),
        statement=(
            "The source study's own classifier is not reproduced here. Its published text and "
            "supplement disagree with each other: a three-feature AUC appears as both 0.74 and 0.76, "
            "and two precision-recall figures, 0.83 and 0.96, are swapped between the two accounts. "
            "The methods mix relative abundance with a read-count adjustment, the cross-validation "
            "is not documented as grouped by household, and the declared code repository was "
            "unreachable. Picking whichever number is more favourable would be the wrong repair."
        ),
        promotion_requirement=(
            "A reproducible locked implementation with source labels, oversampling confined to inner "
            "training folds, original untouched test prevalence retained, and household grouping. The "
            "separate directional index above remains implementable and is what is shown."
        ),
    ),
)

TASKS_BY_ID: Final[dict[str, UnavailableTask]] = {t.task_id: t for t in UNAVAILABLE_TASKS}


# --------------------------------------------------------------------------- #
# invariants
# --------------------------------------------------------------------------- #


def assert_seed_invariants() -> None:
    """Compiler check on this release's inventory (§3.1), run at import.

    A count is a weak check, but a wrong direction or a duplicated slot is
    exactly the sort of error that reads as a plausible result, so it is
    worth catching before anything renders.
    """
    if len(SEED_PANELS) != 9:
        raise DirIndexError(f"v4.2 seed holds nine numeric panels, got {len(SEED_PANELS)}")
    slots = sum(len(p.features) for p in SEED_PANELS)
    if slots != SEED_SLOT_COUNT:
        raise DirIndexError(f"v4.2 seed holds {SEED_SLOT_COUNT} feature slots, got {slots}")
    ids = [p.panel_id for p in PANELS]
    if len(ids) != len(set(ids)):
        raise DirIndexError("panel IDs must be unique")
    for panel in PANELS:
        keys = [f.source_id for f in panel.features]
        if len(keys) != len(set(keys)):
            raise DirIndexError(f"{panel.panel_id}: feature IDs must be unique within a panel")
        measured = [f.key for f in panel.features]
        if len(measured) != len(set(measured)):
            raise DirIndexError(
                f"{panel.panel_id}: two slots resolve to the same measurement, which would count "
                "one observation twice"
            )
        for f in panel.features:
            if f.direction not in (-1, 1):
                raise DirIndexError(f"{panel.panel_id}/{f.source_id}: direction must be -1 or +1")
        if panel.family not in FAMILIES:
            raise DirIndexError(f"{panel.panel_id}: unknown family {panel.family!r}")
    # Alzheimer's AD and dermographism SD must keep their identifiers (§0.2).
    for reserved in ("AD", "SD"):
        clashes = [p.panel_id for p in PANELS if p.panel_id.startswith(f"{reserved}_")]
        if clashes:
            raise DirIndexError(
                f"{reserved!r} is already taken by an existing family; skin panels use ATD/SEBD "
                f"(offending: {clashes})"
            )
    for task in UNAVAILABLE_TASKS:
        if task.family not in FAMILIES:
            raise DirIndexError(f"{task.task_id}: unknown family {task.family!r}")


assert_seed_invariants()


# --------------------------------------------------------------------------- #
# non-scored functional and mechanistic observations (§3.2)
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class Observation:
    """Something the sources measured that carries no disease-score weight."""

    source_id: str
    feature: str
    source_direction: str
    permitted_interpretation: str

    def to_json(self) -> dict[str, Any]:
        return {
            "source_id": self.source_id, "feature": self.feature,
            "source_direction": self.source_direction,
            "permitted_interpretation": self.permitted_interpretation,
            "disease_score_weight": 0.0,
        }


NON_SCORED_OBSERVATIONS: Final[tuple[Observation, ...]] = (
    Observation(
        "E01", "PWY-5005 (biotin biosynthesis)", "increased in discovery, failed to validate",
        "Measured biotin-biosynthesis DNA potential. It failed validation in its own source, so it "
        "carries no disease-score weight and is not a replicated finding.",
    ),
    Observation(
        "E14", "P163-PWY", "decreased",
        "One reported pathway-capacity observation. DNA capacity for a pathway is not a measurement "
        "of any metabolite in stool or blood.",
    ),
    Observation(
        "E14", "PWY-5304", "decreased",
        "One reported pathway-capacity observation; no metabolite level is inferred from it.",
    ),
    Observation(
        "E14", "PWY-7200", "increased",
        "One reported pathway-capacity observation; no metabolite level is inferred from it.",
    ),
    Observation(
        "E04", "KEGG 00500, 02030, 02040", "decreased in all psoriasis versus healthy controls",
        "Source-associated functional potential. Not skin cytokines, and not a measurement of gut "
        "motility or any clinical feature.",
    ),
    Observation(
        "E15", "CAZy GH43 increased; CE4 and CBM50 decreased", "as stated",
        "Source-specific gene-family potential from the same samples as the nine gene families "
        "scored above, so it is the same evidence group rather than an independent replication.",
    ),
    Observation(
        "E01", "vBin_422 and three published contig names", "annotation only",
        "A lead with a name and no sequence. Without the actual assembly and a validated mapping and "
        "reference bundle there is nothing to score, and no reference sequence is generated to fill "
        "the gap.",
    ),
    Observation(
        "E03", "enzyme annotations labelled virulence or resistance",
        "reported by the source as annotation categories",
        "Annotation categories, not validated resistance or virulence gene calls. The source's "
        "clinical framing of them is not imported.",
    ),
)


# --------------------------------------------------------------------------- #
# scoring
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class SkinResult:
    """One panel's outcome, or the typed reason there is not one."""

    panel: Panel
    result: PanelResult | None
    status: str
    reason_codes: tuple[str, ...] = ()
    notes: tuple[str, ...] = field(default=())

    @property
    def computed(self) -> bool:
        return self.result is not None and self.result.computed

    def to_json(self) -> dict[str, Any]:
        if self.result is not None:
            payload = self.result.to_json()
            payload["reason_codes"] = list(dict.fromkeys(list(payload.get("reason_codes", [])) + list(self.reason_codes)))
            payload["notes"] = list(dict.fromkeys(list(payload.get("notes", [])) + list(self.notes)))
            payload["score_revision"] = "directional_fraction_v42"
            return payload
        return {
            "panel_id": self.panel.panel_id,
            "panel_version": self.panel.revision,
            "score_revision": "directional_fraction_v42",
            "profile_family": self.panel.family,
            "label": self.panel.label,
            "contrast": self.panel.contrast,
            "measurement": self.panel.measurement,
            "status": self.status,
            "index": None,
            "index_label": "Research directional pattern index",
            "usable_count": 0,
            "panel_count": len(self.panel.features),
            "coverage": 0.0,
            "interval_coverage": 0.0,
            "feature_mask": [],
            "full_panel_bounds": [0.0, 100.0],
            "bounds_type": "coverage_and_censoring",
            "index_scope": "available_point_feature_mask",
            "bounds_scope": "full_declared_panel",
            "control_percentile": None,
            "case_percentile": None,
            "disease_probability": None,
            "clinical_classification": None,
            "assay_transport": self.panel.assay_transport,
            "source_population_transport": "not_established",
            "disease_specificity": self.panel.disease_specificity,
            "evidence_tier": self.panel.evidence_tier,
            "source_ids": list(self.panel.all_source_ids),
            "counterevidence_ids": list(self.panel.counterevidence_ids),
            "reason_codes": list(self.reason_codes),
            "notes": list(self.notes),
            "feature_results": [],
            "domain_assessment": "limited_metadata_based",
        }


def eligibility(
    panel: Panel,
    *,
    age_years: float | None,
    body_site: str,
    usable_read_pairs: int | None,
    counts_are_pairs: bool,
    qc_passed: bool,
    antibiotics_within_90_days: bool | None,
    mode: str,
) -> tuple[str | None, tuple[str, ...]]:
    """Gates in the specification's order (§7.1), keeping every reason.

    Returns the blocking status, or ``None`` when the panel may be scored,
    together with every reason code that applied. The order matters —
    structure before assay before QC before population before reference — but
    the later reasons are still reported, because "we could not measure this"
    and "you are outside the supported age range" are different facts and a
    reader is owed both.
    """
    reasons: list[str] = []
    blocking: str | None = None

    def block(status: str, code: str) -> None:
        nonlocal blocking
        reasons.append(code)
        if blocking is None:
            blocking = status

    if body_site != "stool":
        block("assay_incompatible", f"body_site_{body_site}_is_not_stool")
    if not qc_passed:
        block("qc_failed", "sample_qc_failed")
    if usable_read_pairs is None:
        block("qc_failed", "usable_read_pair_count_unknown")
    else:
        if not counts_are_pairs:
            block("qc_failed", "read_count_is_reads_not_pairs_so_the_pair_floor_is_unsatisfiable")
        elif usable_read_pairs < MIN_USABLE_READ_PAIRS:
            block("qc_failed", f"usable_read_pairs_{usable_read_pairs}_below_{MIN_USABLE_READ_PAIRS}")
    if age_years is None:
        if mode == "participant":
            block("population_unsupported", "age_unknown_and_unknown_is_not_assumed_eligible")
        else:
            reasons.append("age_unknown_scored_in_research_mode_without_a_population_guarantee")
    elif not (panel.min_age_years <= age_years <= panel.max_age_years):
        if age_years < 18.0:
            block("population_unsupported",
                  f"age_{age_years:g}_is_paediatric_and_adult_disease_numbers_are_not_issued")
        else:
            block("population_unsupported",
                  f"age_{age_years:g}_outside_supported_{panel.min_age_years:g}_{panel.max_age_years:g}")
    if antibiotics_within_90_days is True and mode == "participant":
        block("population_unsupported", f"antibiotics_within_{ANTIBIOTIC_EXCLUSION_DAYS}_days")
    elif antibiotics_within_90_days is None:
        reasons.append("antibiotic_exposure_unknown_and_unknown_is_not_none")
    return blocking, tuple(dict.fromkeys(reasons))


def ko_fractions(measured: Mapping[str, float]) -> dict[str, float]:
    """Normalise gene-family abundances over their own annotation space (§4.1).

    The denominator is the sum across *every* measured microbial gene family
    in the same pinned annotation space, not just the nine this panel selects.
    Normalising within the selected nine would reclose them to 1 and turn any
    change in one into an apparent opposite change in the others.

    Raises when the denominator is not usable, because a zero or negative
    total means the upstream quantity is not what this kernel expects, and
    guessing would be worse than refusing.
    """
    if not measured:
        raise DirIndexError("no gene-family abundances were supplied")
    for key, value in measured.items():
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise DirIndexError(f"gene-family abundance for {key} is not a real number: {value!r}")
        if value != value or value in (float("inf"), float("-inf")):
            raise DirIndexError(f"gene-family abundance for {key} is not finite")
        if value < 0:
            raise DirIndexError(f"gene-family abundance for {key} is negative")
    total = float(sum(float(v) for v in measured.values()))
    if total <= 0.0:
        raise DirIndexError(
            "the gene-family annotation space sums to zero, so no fraction can be formed and every "
            "gene-family slot is unavailable"
        )
    return {k: float(v) / total for k, v in measured.items()}


def _quarantine_faecalibacterium(panel: Panel, catalogue: Sequence[str]) -> tuple[str, ...]:
    """§4.2 case 7: refuse to score an umbrella species beside its own child.

    Historical *F. prausnitzii* was split, and *F. duncaniae* came out of it.
    If the production catalogue does not define the two disjointly then one
    organism would be counted twice with two different directions, so both
    slots are withheld.
    """
    ids = {f.source_id for f in panel.features}
    pair = {"Faecalibacterium_prausnitzii", "Faecalibacterium_duncaniae"}
    if not pair <= ids:
        return ()
    present = {name for name in catalogue if name in pair}
    if len(present) < 2:
        return tuple(sorted(pair))
    return ()


def score_skin_panel(
    panel: Panel,
    abundances: Mapping[str, float],
    reference: FrozenReference | None,
    *,
    catalogue: Sequence[str] = (),
    quarantined: Sequence[str] = (),
    **kwargs: Any,
) -> PanelResult:
    """Score one skin panel, withholding any quarantined slot."""
    withheld = set(quarantined) | set(_quarantine_faecalibacterium(panel, catalogue))
    if withheld:
        trimmed = dict(abundances)
        for key in withheld:
            trimmed.pop(key, None)
            feature = next((f for f in panel.features if f.source_id == key), None)
            if feature is not None:
                for name in (feature.key, *feature.synonyms):
                    trimmed.pop(name, None)
        abundances = trimmed
    return score_panel(panel, abundances, reference, **kwargs)


def skin_indices(
    species_percent: Mapping[str, float],
    cohort: Any | None,
    *,
    age_years: float | None = None,
    body_site: str = "stool",
    usable_read_pairs: int | None = None,
    counts_are_pairs: bool = True,
    qc_passed: bool = True,
    antibiotics_within_90_days: bool | None = None,
    mode: str = "research",
    ko_abundances: Mapping[str, float] | None = None,
) -> dict[str, Any]:
    """Every v4.2 panel for one sample, plus the tasks that produce no number.

    ``species_percent`` is the production taxonomic table in percent, keyed by
    catalogue name; it is divided by 100 exactly once inside the kernel. KO and
    pathway panels are scored only from ``ko_abundances``, never inferred from
    taxa: an organism's presence is not proof of a gene.
    """
    fractions = {k: float(v) / 100.0 for k, v in species_percent.items()}
    catalogue = list(getattr(cohort, "taxa", []) or [])
    matrix = getattr(cohort, "abundance", None)
    out: dict[str, Any] = {}

    for panel in PANELS:
        blocking, reasons = eligibility(
            panel, age_years=age_years, body_site=body_site,
            usable_read_pairs=usable_read_pairs, counts_are_pairs=counts_are_pairs,
            qc_passed=qc_passed, antibiotics_within_90_days=antibiotics_within_90_days,
            mode=mode,
        )
        if blocking is not None:
            out[panel.panel_id] = SkinResult(panel, None, blocking, reasons).to_json()
            continue
        if panel.measurement == "shotgun_ko_fraction":
            # A gene family is measured or it is not. The taxonomic table
            # cannot stand in for it: carrying an organism is not evidence of
            # carrying a gene.
            if not ko_abundances:
                out[panel.panel_id] = SkinResult(
                    panel, None, "features_unavailable",
                    reasons + (
                        "no_pinned_gene_family_annotation_was_supplied",
                        "no_taxa_to_function_imputation_is_permitted",
                    ),
                    notes=(
                        "This panel needs directly annotated gene-family abundances from a pinned "
                        "annotation pipeline, normalised across every measured microbial gene family "
                        "in that same space. The taxonomic table cannot stand in for it.",
                    ),
                ).to_json()
                continue
            try:
                ko_fractions(ko_abundances)
            except DirIndexError as exc:
                out[panel.panel_id] = SkinResult(
                    panel, None, "assay_incompatible",
                    reasons + ("gene_family_denominator_unusable",),
                    notes=(str(exc),),
                ).to_json()
                continue
            out[panel.panel_id] = SkinResult(
                panel, None, "reference_unavailable",
                reasons + ("no_compatible_frozen_gene_family_reference",),
                notes=(
                    "Gene-family abundances were measured and normalised, but no frozen reference "
                    "distribution exists for them in a compatible annotation space. The source's own "
                    "ten controls do not meet this build's minimum of twenty per feature.",
                ),
            ).to_json()
            continue
        if not catalogue or matrix is None:
            out[panel.panel_id] = SkinResult(
                panel, None, "reference_unavailable", reasons + ("no_compatible_frozen_reference",),
            ).to_json()
            continue
        frozen = freeze_reference(panel, catalogue=catalogue, matrix=matrix, bundle_id="refs/taxonomic_cohort")
        result = score_skin_panel(panel, fractions, frozen, catalogue=catalogue)
        out[panel.panel_id] = SkinResult(
            panel, result, result.status, reasons,
        ).to_json()

    out["unavailable_tasks"] = {t.task_id: t.to_json() for t in UNAVAILABLE_TASKS}
    out["non_scored_observations"] = [o.to_json() for o in NON_SCORED_OBSERVATIONS]
    out["families"] = {k: v.to_json() for k, v in FAMILIES.items()}
    out["release"] = RELEASE
    out["no_family_fusion"] = (
        "Panels within one family are reported side by side and never averaged into a family score. "
        "Several panels do not make several diseases, and one feature appearing in several panels is "
        "one observation."
    )
    return out
