# OpenBiota — BUILD_SPEC_v0.8.3

**Release objective:** extend the existing OpenBiota report with a substantially broader set of interpretable measurements, substrate-specific functions, ecological context, evidence-linked actions, and longitudinal reporting.

**Baseline:** the running v0.8.2 application and its existing results contract. **Research cutoff:** 22 September 2026. **Change policy: additive; broaden and correct new action retrieval/ranking as specified in sections 8, 11 and 19.** This document is the complete implementation handoff; its research, source links, data requirements, algorithms, feature inventory, and acceptance criteria are embedded here.

## 1. Binding scope and preservation contract

Implement additions within the existing report. Preserve its calculations, numerical results, scientific interpretation and capabilities. The presentation changes explicitly specified here are authorized: integrate new readings into existing sections, broaden the functional section title, replace the lower part of page 3 with contents/navigation, and consolidate inputs and assumptions in report section 25.

Preserve existing disease scores, GMWI2, donor matching, strain analysis, pathogen screening, biofilm scores, mycobiome scores, functional panels, and treatment-evidence sections. New modules must not feed into those existing calculations implicitly.

At implementation start, capture the current production commit, configuration, databases, model artifacts, and representative `results.json` files. Create a protected-output manifest from the actual repository, including every existing metric ID. Run identical inputs through the baseline and the updated application. Existing metric values, status labels, component values and calculation fingerprints must remain unchanged for identical resolved inputs. Retain scientific explanations while allowing the heading, grouping, relocation and navigation changes specified in section 12 of this specification. Move existing explanatory content to its designated destination rather than duplicating it or losing it. Regenerate page numbers, section numbers, contents links and bookmarks after layout.

Integrate new work into the existing application using its established modules, result schema, registries, report templates and cache conventions. Add backward-compatible fields and capabilities where needed; do not rewrite existing cache objects. If a proposed addition requires changing a protected algorithm, implement it as a separate, clearly named companion measurement instead. An unresolved issue in an existing module is not permission to redesign that module under this specification.

### 1.1 Baseline capabilities and input verification

Begin with a repository and installed-asset inventory for the running OpenBiota application. Identify the existing organism, strain, functional, disease, pathogen, biofilm, mycobiome and intervention modules, their canonical result objects, and their data dependencies. Reuse working capabilities and extend their existing interfaces.

The baseline configuration includes MetaPhlAn and a **GTDB r232 / Sylph 1.0.0** detection lane. Verify the actual installed versions, reference manifests and execution states before implementing adapters or additional analyses. Capture representative structured results and rendered OpenBiota reports as preservation fixtures.

Verify specimen IDs, collection dates, FASTQ hashes, library preparation and pipeline fingerprints before combining data or calculating longitudinal changes. Acceptance is determined by the feature contracts, primary research, versioned reference assets and tests defined in this document.

### 1.2 What counts as completion

Every feature below must have one of four explicit implementations:

1. **Existing capability, newly surfaced:** read the existing measurement and add a useful view without changing it.
2. **New sequence-derived measurement:** implement the documented analysis, quantify its evidence, and render its result.
3. **New research model or mechanism panel:** implement a transparent experimental output with its assumptions and supporting components visible.
4. **External-assay integration:** accept and display a measured laboratory result when stool DNA cannot supply that measurement.

An empty card, a literature paragraph without its specified computation, or a `TODO` source URL does not complete a sequence-derived feature. Conversely, inventing a host hormone concentration or an RNA-virus negative result from a DNA-only assay does not complete an unavailable measurement. Implement the capability and a precise data-availability state, so a user can see what additional input would unlock it.

## 2. Additive feature map

| ID | Addition | Existing capability to retain | Required new output |
|---|---|---|---|
| A01 | Fibre and carbohydrate utilization | Existing SCFA panels and fibre advice | Separate cellulose, general starch, resistant starch, chitin, pectin, inulin, FOS, GOS, XOS, IMO, beta-glucan, arabinoxylan, guar/galactomannan, and lactose capacities |
| A02 | Fermentation network | Existing butyrate and propionate | Acetate; lactate production and consumption; succinate production and consumption; hydrogen production and consumption; cross-feeding links |
| A03 | Complete vitamin pathway overview | Existing B2, B7, B9, B12 and K2 | B1, B3, B5 and B6; synthesis versus salvage; one complete nine-vitamin dashboard |
| A04 | Protein and nitrogen metabolism | Existing BCAA and urease readings | Proteolysis, amino-acid fermentation, non-urease ammonia routes, branched-chain fatty-acid capacity |
| A05 | Expanded neuroactive metabolism | Existing GABA, histamine, tryptamine and tyrosine-decarboxylase results | GABA breakdown; reaction-specific catecholamine handling; experimentally supported serotonin routes; acetylcholine and norepinephrine evidence/assay cards |
| A06 | Polyphenol and plant-compound conversion | Existing polyphenol organism group and readings | Separate urolithin reactions, S-equol conversion, glucosinolate conversion, and research links to dietary substrates |
| A07 | GLP-1-related microbial mechanisms | Existing SCFA and bile-acid readings | A component-based mechanism panel, including known microbial proteins where sequence support exists; measured GLP-1 import |
| A08 | Mucin, LPS and bile-acid detail | Existing mucus, endotoxin and bile-acid panels | Substrate-specific mucin enzymes, lipid-A modification evidence, and additional bile transformations without altering existing scores |
| A09 | Ecology companion dashboard | Existing Shannon, evenness, richness and GMWI2 | Gut community type/context, dominance, Hill diversity, composition ratios, aerotolerance balance, functional redundancy and longitudinal resilience |
| A10 | Complete organism explorer | Existing inventory, interpretation and strain lanes | Searchable lineage/alias views, all named screening targets, contributor lists, confidence/coverage, probiotic-strain distinctions, full printable inventory |
| A11 | Resistance overview | Existing resistance detections | Class-level abundance, determinant richness, mechanism and carrier-evidence tables |
| A12 | Consolidated action planner | Existing treatment evidence and recommendations | Deduplicated food, fibre, probiotic, supplement, herb and lifestyle options tied to every relevant existing/new finding, with broad experimental evidence, exact identities and delivery/tradeoff explanations |
| A13 | Longitudinal and context views | All existing results | Last-three-results charts, comparable-run changes, symptom/context tiles, exposure timeline and optional intake data |
| A14 | Your information, inputs and assumptions | Existing metadata inputs, defaults and sequence analyses | One combined report section 25 covering supplied facts, effective assumptions, optional questionnaire inputs and useful laboratory imports linked to relevant findings |
| A15 | Integrated layout and navigation | Existing section structure and readings | New readings in their existing domain sections; clickable major-section contents on page 3; linked details, complete metric index, evidence links and printable actions |
| A16 | Personalized synbiotic options | Existing metabolism, action and matching algorithms | Required strain–substrate evidence matching and transparent candidate ranking; optional reproducible metabolic scenarios with explicit coverage and uncertainty |

The detailed inventories later in this document are normative. Do not implement only the headline rows. A01–A16 identify implementation capabilities, not sixteen new report sections. The report placement matrix in section 12 is binding.

## 3. Integration architecture

### 3.1 Repository integration

Inspect the actual repository and extend its existing CLI, result schemas, report templates, functional analysis modules, taxonomy adapters, reference-cohort builder, intervention registry and tests. Map the A01–A16 capabilities onto those established components. Add a module only when a capability needs one, placing it alongside related code and following existing naming conventions.

The release version identifies this specification and relevant data/model versions; it does not prescribe a source-code directory or a new result namespace. Do not create a separate version-specific directory tree or duplicate the application pipeline. Use existing interfaces and shared services, keeping new calculations independently testable and preserving existing behavior. Concrete file paths must come from the repository inspected by the coding agent.

The production handoff is this single document. The coding agent should create normal repository code, registries, tests, downloaded data, and cache manifests during implementation; those are application assets, not missing attachments to this specification.

### 3.2 Inputs and incremental execution

Prefer validated existing data in this order: structured result objects; versioned intermediate tables; alignment/gene-count caches; retained host-filtering logs; FASTQs. PDFs are display artifacts and must not become the main calculation input when structured files exist.

For each addition, emit a dependency manifest containing input hashes, database release, sequence/reference hashes, parameter hash, program version, feature-manifest version, and reference-cohort fingerprint. Cache reuse requires all relevant dependencies to match. A new CAZyme panel does not invalidate unrelated existing caches.

Support three execution modes:

- `cached`: create every possible new view from existing structured data; explicitly list additions needing further analysis.
- `incremental`: run only missing sequence analyses, retaining existing outputs unchanged.
- `full-extension`: calculate the complete extension from eligible reads and cached baseline results; still do not rerun or recalibrate protected modules unless their normal, unchanged workflow already requires it.

Missing FASTQs must not prevent the existing report from rendering. It must produce the available additions and an exact list of unavailable analyses. Production completion, however, requires an integration fixture with FASTQs that exercises all implemented sequence-analysis modules.

### 3.3 New result contract

Every new measurement has a stable ID and the following fields. Fields may be extended, never overloaded:

```json
{
  "schema_version": "openbiota.extension-metric/1.0",
  "metric_id": "ext083.function.acetate.pta_ackA",
  "label": "Acetate formation: phosphotransacetylase/acetate kinase route",
  "kind": "genetic_capacity",
  "state": "measured",
  "value": 12.4,
  "unit": "fragment_RPKM",
  "denominator": "quality_passed_nonhost_fragments",
  "direction": "context_dependent",
  "reference_percentile": null,
  "reference_id": null,
  "reference_state": "not_available",
  "score_0_100": null,
  "score_definition_id": null,
  "analytical_confidence": "supported",
  "evidence_maturity": "biochemical_route",
  "pathway_completeness": 1.0,
  "completeness_scope": "community_assembled",
  "linkage_evidence": [],
  "assessable_fraction": 1.0,
  "components": [],
  "contributing_taxa": [],
  "source_ids": [],
  "method_id": "ext083.route_capacity/1.0",
  "input_fingerprint": "sha256:example",
  "limitations": ["Genetic capacity does not measure acetate concentration."],
  "action_links": [],
  "legacy_metric_links": []
}
```

The numerical example is a schema illustration, not a result for any supplied sample.

`kind` must distinguish `taxon_abundance`, `gene_abundance`, `genetic_capacity`, `community_composition`, `reference_percentile`, `experimental_index`, `model_prediction`, `external_lab_measurement`, and `evidence_context`.

`state` must distinguish `measured`, `partial`, `not_detected_above_assay_threshold`, `insufficient_coverage`, `not_assayed`, `unsupported_by_assay`, `requires_external_result`, and `not_applicable`. Zero, null, and not assessed are different states. A reference percentile is not a pathway-completeness percentage, disease probability, concentration, or probability of treatment benefit.

`direction` is one of `higher_favourable_in_context`, `lower_favourable_in_context`, `interval_favourable_in_context`, `context_dependent`, or `descriptive`. It must cite the context that supports the direction. New generic abundance measurements default to descriptive, not automatically green when high.

Every `contributing_taxa` row must state the evidence of assignment: direct strain-specific sequence, uniquely mapped gene, supported contig linkage, taxonomically classified MAG, ambiguous reference hit, or organism-level association. Co-detection of a species and a gene is not sufficient to claim that species carries that gene in the sample.

### 3.4 New gene and reaction registry

A gene symbol is not a unique biochemical identity. Each new marker record must include sequence accession and version, amino-acid sequence hash, source publication, organism/strain, assay evidence, substrate, product, reaction direction, orthology family, diagnostic residues when needed, negative homologs, and qualifying context.

Use accession/reaction IDs as joins. Never bind a new metric by a display-name substring. In particular, the new urolithin, equol, and imidazole-propionate reaction IDs must be distinct; `urdA` is not a urolithin-conversion marker. This is a requirement for the new modules and does not authorize rewriting existing outputs.

For every marker family, create a positive set and a difficult negative/near-homolog set. A reference sequence added merely because its annotation contains a desired word is an unvalidated candidate until specificity has been checked.

## 4. New functional analysis engine

### 4.1 Reference construction and quantification

1. Reuse existing gene alignments only if their reference sequences and parameters actually cover the new markers.
2. Build a separate extension protein database from experimentally characterized proteins, curated pathway definitions, and reference-genome proteins annotated with the selected tools.
3. For carbohydrate enzymes, run dbCAN/HMM and substrate/PUL annotation on **full-length reference or assembled proteins**, then map reads to the annotated protein catalogue. Do not run a full-length protein-domain workflow directly on raw FASTQs and treat a short conserved hit as substrate-specific proof.
4. Use competitive alignment against both positive and homologous-negative sequences. Count paired reads as one fragment when they originate from the same molecule. Preserve alternative mappings and assign fractional support rather than counting the same fragment in full against every homolog.
5. Report read support, covered residues/regions, sequence identity, alignment length, reference specificity, and pathway step coverage. A short hit to a ubiquitous enzyme may support a broad family while remaining insufficient for the specific reaction.
6. Optional assembly/MAG analysis adds carrier and pathway-colocalization evidence. Absence from a fragmented assembly does not override a supported read-level detection.

For a reference gene of effective length `L_g` nucleotides and assigned fragments `n_g`, define the new lane's abundance as:

```text
fragment_RPKM(g) = 1e9 * n_g / (L_g * N_nonhost_fragments)
```

Store actual effective lengths and multi-mapping treatment. A protein reference uses its documented CDS length or a defined nucleotide-equivalent length (`3 × amino-acid length`, with stop-codon treatment fixed); never use amino-acid length with a nucleotide-length unit. Partial assembled genes require an explicit effective-length model and must not inflate abundance simply because the contig is short. If reporting genome-equivalent copies, use a documented multi-single-copy-marker normalization model calibrated for the domains included, with a separate unit and denominator. Do not silently replace the legacy functional normalization. Reads per kilobase, fragments per kilobase, bacterial genome equivalents, and fungal genome equivalents are not interchangeable.

For an AND step requiring multiple subunits, use the minimum abundance, applying subunit stoichiometry only when established. For an OR step, sum non-overlapping alternative reactions, with ambiguous evidence counted once. The minimum across required steps is an **engineering bottleneck-abundance proxy**, not measured flux or a validated biochemical production rate. Publish the step values and completeness so an isolated marker does not masquerade as a complete route. `completeness_scope` distinguishes `genome_linked`, `community_assembled`, and `reference_inferred`; record linkage evidence. Community-level cross-feeding is distinct from all steps coexisting in one linked genome.

### 4.2 New reference percentiles and experimental scores

New percentiles require reference samples processed using the exact new method. Existing percentiles remain unchanged. A sample lacking an eligible reference still receives supported raw capacities and evidence, rather than losing the entire feature.

Use participant-level, study-aware reference sets with phenotype and collection metadata. Do not label every public stool sample healthy. Store exclusions and retain age, geography, sex where available, medications, assay, and processing information. Reference cohorts restricted to a particular disease study may supply a named research comparison but are not automatically a general healthy range.

For a reference distribution `R` and value `x`, use a documented midrank empirical percentile:

```text
P(x; R) = 100 * (count(R < x) + 0.5 * count(R == x)) / len(R)
```

Use one participant vote or fixed study-balanced weights; do not allow repeated visits or a single large study to silently dominate. Store the choice. Quantile bands are descriptive unless an outcome-linked normal interval has been established. For sparse genes, show prevalence and a positive-carrier percentile separately, so many zero values do not create a misleading low/high ranking.

A normalized 0–100 experimental index must disclose its component formula, coverage, and version. It may be useful without being a clinical test. Never fabricate a percentile without a reference, hide missing components by replacing them with healthy values, or change existing scores to accommodate the extension.

### 4.3 Functional display contract

Every new functional card includes: what the microbes can potentially do; which route was measured; raw result and unit; coverage; reference context if available; the main supported contributors; relevant foods or substrates; actionable evidence; and what would directly measure the endpoint. Sequence-derived capacity and an imported metabolite concentration can appear beside each other but remain separate measurements.

## 5. Required functional additions and their implementation

Source IDs F01–F43 resolve to the embedded functional source register below. Implement these modules independently of the existing 25 panels. Where a reaction is already measured, expose its cached evidence through a new view rather than running a duplicate analysis.

### 5.1 Dietary substrate panel — A01

Build these **14 distinct substrate views**, with shared evidence explicitly linked:

| New ID suffix | Substrate | Required discrimination | Useful food context |
|---|---|---|---|
| `carb.cellulose` | Cellulose | Substrate-validated cellulases/accessory enzymes; distinguish other beta-glucans | Vegetables, bran and intact plant material |
| `carb.starch_general` | General starch | Amylolytic enzymes and uptake; broad comparator for the next row | Starchy foods; cooking/processing context |
| `carb.resistant_starch` | Resistant starch | Characterized resistant-starch utilization systems; record supported physical type | Legumes, unripe banana, cooked/cooled starches, named RS preparations |
| `carb.chitin` | Chitin | Chitinase plus relevant downstream utilization; distinguish peptidoglycan | Mushrooms and other chitin-containing foods, with dietary context |
| `carb.pectin` | Pectin | Pectin-specific lyases/hydrolases, esterases and relevant PULs | Apples, citrus, vegetables |
| `carb.inulin` | Inulin | Long-chain fructan substrate evidence and uptake context | Chicory, alliums, artichoke-related foods |
| `carb.fos` | Fructooligosaccharides | Short-chain fructan evidence; share but do not duplicate inulin genes | FOS-containing foods and defined preparations |
| `carb.gos` | Galactooligosaccharides | GOS linkage/substrate evidence and transport, beyond generic beta-galactosidase | Defined GOS preparations; distinguish them from all legume carbohydrates |
| `carb.xos` | Xylooligosaccharides | Xylosidases and transport/chain-processing context | Defined XOS; relevant whole-grain substrates |
| `carb.imo` | Isomaltooligosaccharides | Alpha-glucoside linkage and chain-length specificity | Defined IMO preparations; formulation matters |
| `carb.lactose` | Lactose | Characterized lactose hydrolysis/transport routes | Dairy lactose and fermented-food context |
| `carb.beta_glucan` | Beta-glucans | Linkage-specific enzymes; cereal and fungal beta-glucans remain distinguishable | Oats/barley versus fungal preparations |
| `carb.arabinoxylan` | Arabinoxylan | Xylan backbone plus arabinose side-chain processing | Wheat/rye bran and specified extracts |
| `carb.galactomannan` | Galactomannan/guar | Mannan backbone, side-chain processing and characterized substrate evidence | Guar and partially hydrolysed guar gum |

Use dbCAN family, subfamily, substrate and PUL resources [F01–F05]. Required assets include `CAZy.dmnd`, `dbCAN.hmm`, `dbCAN-sub.hmm`, `fam-substrate-mapping.tsv`, `PUL.dmnd`, and the dbCAN-PUL substrate metadata. On assemblies, retain CGC/PUL synteny, accessory enzymes, transporters, and secretion/localization evidence where available. On low-depth samples, quantify curated reference genes directly and mark unavailable genomic context separately.

The reviewed run_dbcan V5 code snapshot is `c23f0d08d7e2678feca7b485bed4b3a4ce7b694a`. The official server lists HMMdb v15 dated 6 August 2026, based on CAZy 24 July 2026. Fetch a **mutually compatible complete snapshot**, not a mixture of this HMM with an older substrate mapping. The AWS `db_v5-2_9-13-2025` collection is a different dated snapshot. Record the actual selected release and checksums. Use the version's documented `run_dbcan database --db_dir ...` and substrate-analysis interfaces after testing its supplied examples.

Specificity fixtures are mandatory: GH13 alone supports general amylolysis, not every resistant starch; GH32 alone does not separate FOS from inulin; GH2/GH42 alone do not establish specific GOS utilization; GH18 alone does not demonstrate an intact chitin-utilization pathway. A shared gene can appear as linked evidence in two substrate cards, but its fragments cannot be counted twice in an aggregate.

A new **fibre opportunity table** joins substrate capacity, likely downstream fermentation routes, food examples, reported tolerance, and evidence. Low capacity does not prove a food cannot be tolerated; high capacity does not prove it will increase butyrate or improve symptoms. Show the relevant uncertainty and offer alternative substrates. Microbial lactose utilization does not measure the person's small-intestinal lactase.

### 5.2 Fermentation and cross-feeding — A02

Add the following explicit measurements and a compact substrate → intermediate → terminal-product diagram:

| New metric | Sequence/route definition | Interpretation |
|---|---|---|
| Acetate formation | `pta` + `ackA` route and other experimentally characterized alternatives; separately reconstruct acetogenic Wood–Ljungdahl metabolism | Potential acetate formation; ordinary AMP-forming `acs` supports acetate assimilation, not formation by itself; distinguish ADP-forming enzymes and Wood–Ljungdahl CODH/ACS |
| Lactate formation | Curated D- and L-lactate-forming enzymes with reaction context | Report stereoisomers where resolved |
| Lactate utilization | Validated lactate-utilizing routes, including relevant acrylate/propionate and butyrate cross-feeding systems | A complementary sink, not a measured subtraction from production |
| Succinate formation | Supported fermentative route modules | Potential intermediate supply |
| Succinate utilization | Validated succinate-to-propionate and other specific consumption routes | Potential intermediate consumption |
| Hydrogen production | Hydrogenase classes, particularly direction-resolved [FeFe] systems | No universal adverse polarity |
| Hydrogen consumption | Hydrogenotrophic methanogenesis, acetogenesis and relevant respiratory systems, retaining their distinct endpoints | Link to existing methane/sulfide results without changing them |

Use route reconstruction and validated sequences from F08–F14 and F24–F27. Preserve the existing butyrate, propionate, methane, H2S and TMA scores; add a separate pathway breakdown. For propionate, distinguish succinate, acrylate and propanediol chemistry. For sulfur, distinguish sulfate/sulfite respiration, taurine/isethionate, sulfoacetate and cysteine routes. The `IslA` enzyme requires its appropriate activase/context; a generic glycyl-radical enzyme is not enough. Do not equate *Bilophila* abundance with every sulfur pathway.

Use HydDB class definitions and the 2025 hydrogenase study [F24–F25]; the study's public datasets are **PRJEB45397** and **PRJEB70412**. Its high hydrogen-producing potential in healthy communities is a reason to show production and consumption separately, not to label all hydrogen capacity bad.

The cross-feeding graph displays **supported biochemical opportunities**. Solid edges require supported reaction evidence; dashed edges indicate literature-compatible hypotheses. It must not present an unmeasured flow rate, assume substrates are available, or claim that a taxon will engraft because it fills one node. The existing donor matcher remains unchanged.

### 5.3 Nine-vitamin dashboard — A03

Add B1, B3, B5 and B6. Reuse B2, B7, B9, B12 and K2 unchanged. The new dashboard includes synthesis, salvage and uptake columns; a transport gene does not count as synthesis.

| Vitamin | New reconstruction rules |
|---|---|
| B1 / thiamine | Reconstruct thiazole and hydroxymethylpyrimidine branches plus coupling/activation; distinguish salvage of thiamine and its precursors |
| B3 / niacin-related cofactors | Reconstruct aspartate → quinolinate → NaMN → NAD pathways separately from nicotinamide/nicotinate salvage |
| B5 / pantothenate | Reconstruct pantoate and beta-alanine branches plus pantothenate ligation; downstream CoA utilization is not a second production pathway |
| B6 | Accept documented DXP-dependent PdxA/PdxJ and DXP-independent PdxS/PdxT alternatives; do not require both or score one isolated homolog as a complete route |

Import the curated pathway/gene inventories from the supplements to Magnúsdóttir 2015 and Rodionov 2019 [F06–F07]. Record every imported gene, alternative and source table. Their genome reconstructions supply route definitions and test genomes, not sample-specific evidence of production.

For existing vitamins, add an explanatory route-detail view only. B12-related genes may support different corrinoids or salvage, and K2-related genes do not automatically identify a specific menaquinone chain length. The new dashboard reports microbial capability, not serum vitamin status, host absorption, or a reason to stop a prescribed supplement.

### 5.4 Protein, nitrogen and aromatic metabolism — A04

Add separate modules for:

- Proteolytic/peptidase capacity, with MEROPS family assignments, localization and substrate specificity where established.
- Peptide transport and amino-acid utilization.
- Oxidative and reductive amino-acid fermentation, including supported Stickland routes.
- Branched-chain amino-acid fermentation to branched-chain fatty acids, distinct from existing BCAA biosynthesis.
- Non-urease ammonia-generating reactions and ammonia assimilation. Retain urease as its own existing result.
- Aromatic precursor routes such as phenylacetate, linked to existing indole, IPA and p-cresol readings without replacing their chemical separation.

Use F10–F14, curated reaction resources and experimentally characterized enzymes. Generic peptidase abundance is not a validated measurement of dietary protein digestion. Protein breakdown, BCAA synthesis and BCAA fermentation are three different concepts. Microbial phenylacetate is a precursor; host conjugation produces circulating PAGln. Likewise, microbial TMA formation is distinct from host TMAO production. Provide these relationships as additional explanatory links.

### 5.5 Neuroactive pathways — A05

**GABA breakdown is a required measured addition.** Inspect the existing OpenBiota gene table and reuse supported `gabT` observations with their provenance. Reconstruct the GabT/GabD route and validated alternatives, including the 4-aminobutyrate-to-butyrate route where supported. Display production and degradation side by side. A difference between two percentiles is not net flux. If showing a balance ratio, use compatible raw units, state its formula and show both components.

Use published Gut Brain Module definitions [F15–F16] as a source of reaction logic. Verify the licensing of both definitions and software: the reviewed omixer wrapper and bundled Java program have different terms. An independent scorer over lawfully usable definitions avoids making a noncommercial binary a mandatory production dependency.

Preserve the native AND/OR grammar and identifier namespaces: the released modules mix KO, COG, TIGRFAM and legacy eggNOG identifiers. Relevant definitions include MGB019 GABA degradation, MGB020–022 GABA synthesis alternatives, MGB043–046 acetate synthesis and MGB047 acetate degradation. MGB008/MGB013 historical neurotransmitter labels do not by themselves establish substrate-specific activity. Do not replace unresolved old identifiers with similarly named current proteins.

Add these cards:

| Card | Implementable content | Endpoint boundary |
|---|---|---|
| Dopamine-related microbial metabolism | Link the existing tyrDC reading with its current status; label specific L-DOPA conversion only with substrate-specific evidence, and add validated downstream routes | Not brain dopamine production or a medication-dose instruction |
| Serotonin-related microbial metabolism | Include the 2025 *L. mucosae* / *L. ruminis* consortium result and its 5-HTP-dependent chemistry [F17]; use exact strain/enzyme evidence if verified | Do not substitute tryptophan-to-tryptamine markers or infer circulating serotonin |
| Acetylcholine | Show experimentally supported strain/reaction evidence and compatible imported measurements; implement a sequence score only for a validated specific panel | Generic acetyltransferases are not an acetylcholine health assay |
| Norepinephrine | Separate microbial handling/deconjugation and host-response evidence; display applicable measured external results | No generic bacterial abundance-to-norepinephrine conversion |
| Histamine detail | Link existing production reading; add a specific degradation route only where its enzyme activity and sequence discrimination are established | Potential production/consumption is not a diagnosis of histamine intolerance |

F17 supplies consortium evidence and a traceable candidate decarboxylase, not a validated serotonin-specific gene panel. Implement the consortium/context card and candidate-enzyme evidence with exact identities; use the accession, normalized-sequence hash and unresolved construct/substrate-specificity details in F17. Do not upgrade a tryptamine-producing enzyme to a serotonin detector merely because a homolog is retrievable. F18 supplies contrasting commensal evidence and reinforces strain/condition specificity.

### 5.6 Polyphenols and glucosinolates — A06

Add distinct new metric IDs and independent reference distributions:

1. `polyphenol.urolithin_9_dehydroxylation`: use the experimentally characterized **ucdCFO** complex; **GenBank PQ855390.1** supplies a validated operon. Require reaction-specific evidence for the UcdC, UcdF and UcdO components, with locus coherence reported where assembly permits. This supports the urolithin C → urolithin A step, not the entire ellagitannin-to-urolithin pathway. Include precursor and alternative regioselective routes from F19–F20 as separately identified reactions.
2. `polyphenol.equol_daidzein_conversion`: reconstruct validated DZNR/dzr, DHDR/ddr, THDR/tdr and relevant DDRC chemistry using F21. Aliases such as `eqlA/B/C` require source/strain-specific mapping. General reductase domains are insufficient.
3. `diet.glucosinolate_isothiocyanate_conversion`: seed with the characterized *B. thetaiotaomicron* BT2159–BT2156 locus [F23]. The experimentally supported in-vitro core requires **BT2158 AND (BT2156 OR BT2157)**. Retain the full-locus context and distinguish this bacterial route from plant myrosinase and other bacterial alternatives.

Add food-source links: ellagitannin-containing foods, soy isoflavones, and cruciferous vegetables respectively. Show processing/formulation context and route coverage. DNA evidence does not by itself establish a measured urolithin/equol producer metabotype; substrate-challenge metabolomics can be imported into the same section as a separate measurement.

Use UrdA and near homologs as negative controls for the new urolithin panel [F22]. Its substrate is urocanate and its product is imidazole propionate. Residue-based discrimination must use aligned reference numbering and covered diagnostic positions, not a raw-read position. Existing legacy IDs and outputs are not repurposed.

### 5.7 Mucin, lipid A and bile-acid detail — A08

**Mucin:** add characterized mucin-glycan CAZyme/PUL modules, including relevant sialidase, fucosidase, sulfatase and downstream glycan-processing evidence [F41–F42]. The 2025 resource characterizes 66 enzymes; the July 2026 colonic-mucin study adds substrate-specific sulfatase evidence and has a September correction that must be ingested. Preserve gastric versus colonic substrate context. Use substrate-specific annotations and distinguish mucin-associated enzymes from dietary-fibre enzymes. Link the existing mucin-associated organism group. A physiological mucin user is not automatically damaging the mucus barrier; genetic capacity is not an image of mucosal erosion.

**Lipid A/LPS:** add a structural-evidence table alongside the existing hexa-LPS group. Record supported lipid-A acylation/modification enzymes, species/strain context, and the known versus unresolved predicted structure. Do not infer exact TLR4 stimulation from total Gram-negative abundance, a generic `lpx` hit, or one modification gene. When the structure cannot be resolved, show the gene evidence and uncertainty rather than a fabricated endotoxin concentration.

**Bile acids:** add BSH substrate preference, bai pathway completeness, position/stereospecific HSDH transformations, and microbial reconjugation as distinct components [F33–F36]. BSH can have both hydrolase and acyltransferase activity, and some HSDHs act on conjugated substrates. Build a reaction graph rather than forcing all chemistry through one universal linear sequence. Preserve the existing bile-acid scores.

BileActome offers 27 HMM families and a reviewed code snapshot `2641f51eed97c6880adaedf87ab374033f84850c` [F36]. Use its sequences/models as candidate annotations after validation against human-gut biochemical examples. Its rumen validation population must not supply human healthy percentiles. The authors' E-value rule is an analysis starting point, not a clinical threshold.

### 5.8 Additional useful capacities

Implement these as clearly separated research additions after the required panels above:

| Module | Required chemistry and source | Output |
|---|---|---|
| Anaerobic urate degradation | Curated anaerobic clusters from F28–F29; distinguish oxygen-requiring uricase and purine salvage | Known-route abundance, completeness and contributors; not serum urate |
| Polyamine pathways | Putrescine, agmatine and spermidine; canonical and validated carboxyaminopropylagmatine alternatives [F30–F31] | Separate capacities, with context-dependent interpretation |
| Glutathione synthesis | `gshA` + `gshB` or validated fused `gshF`; distinguish utilization/recycling [F32] | Known-route potential, acknowledging uncharacterized alternatives |
| Glucuronidase substrate detail | Validated loop/domain and catalytic-motif classes [F39] | New substrate-specific table; no automatic estrogen or drug exposure estimate |
| Drug-metabolism evidence | Experimentally characterized microbial reactions, exact drug and enzyme [F40] | Mechanism/context card with the existing report unchanged; no dose adjustment |

Implement these additional capabilities using the same evidence, coverage and reference contracts as the core panels. “Research” is a visible evidence level, not permission to use untraceable markers.

## 6. Ecology, community type and longitudinal additions

### 6.1 Companion diversity and dominance views

Retain all existing Shannon, evenness and richness results. Add companion values calculated from one explicitly declared native abundance vector `p`, closed to sum one over its stated accepted inventory:

```text
H_natural = -sum(p_i * ln(p_i))
effective_Shannon_species = exp(H_natural)
inverse_Simpson = 1 / sum(p_i^2)
Gini_Simpson = 1 - sum(p_i^2)
dominance_top1 = max(p_i)
dominance_top5 = sum(five largest p_i)
```

These make diversity easier to understand. They do not replace an existing entropy calculation with a different log base or different detected-species set. Display detection depth, catalogue, abundance filter, and unresolved fraction. Require compatible assay and processing fingerprints before interpreting numerical changes as biological changes.

Add descriptive ratios from a single compatible abundance vector:

- Proteobacteria/Pseudomonadota : Actinobacteriota.
- Prevotella sensu the declared taxonomy : Bacteroides.
- Fusobacterium : Faecalibacterium.
- Link the existing Firmicutes/Bacillota : Bacteroidota result unchanged.

Show numerator, denominator, taxonomy membership and units. A ratio is dimensionless; if formatted as a percentage it must explicitly mean `100 × numerator/denominator`, not community abundance. Zero denominator yields undefined/censored status. No universal healthy cutoff is assigned merely to fill a green/red gauge, and no action may recommend increasing *Fusobacterium* just to move a ratio.

### 6.2 Gut community type — A09

Add a composition chart and a research community-type classifier. The always-available output is a genus/family composition view using non-overlapping partitions. Build the classifier from the enterotype methods and datasets in [E01–E03], with explicit feature transformations, clustering parameters and frozen model artifacts.

1. Assemble a frozen reference using one abundance method and a fixed genus universe, plus unresolved/other. Normalize each sample to one.
2. Compute Jensen–Shannon divergence with natural logs and cluster using PAM on `sqrt(JSD)`.
3. Evaluate `k = 2..6` with study-stratified bootstrap stability and held-out silhouette. Freeze the selected model, medoids, features and evaluation report. If discrete clustering is unstable, return continuous composition and distances instead of forcing a category.
4. Assign new samples to the nearest frozen medoid. Report all distances. Optional soft similarities are `softmax(-distance/tau)`, with `tau` selected on held-out reference data. Label them similarities, not disease probabilities.
5. Name groups using their actual distinguishing taxa. Close or unstable assignments are mixed/intermediate. Repeated-sample and read-bootstrap support may be shown where those data exist.

A Bacteroides-dominant result does not prove excessive protein/fat intake. A relative configuration resembling a published low-load enterotype does not establish low microbial biomass without an absolute-load assay. Keep the composition useful without adding unsupported metabolic claims.

### 6.3 Aerotolerance and oxygen-related ecology

Add `aerotolerant_fraction`, `MAPI`, and trait coverage [E04–E07]. Let `A` be supported aerotolerant abundance, `N` strict-anaerobe abundance, and `U` unresolved abundance, all from one lane:

```text
aerotolerant_fraction = A / (A + N)
trait_coverage = (A + N) / (A + N + U)
MAPI = ln(A / N)              # only when A > 0 and N > 0
```

If `A+N=0`, return no supported traits. With one component zero, the fraction remains computable, but MAPI is censored; do not insert an arbitrary pseudocount. Display the contributor table and the unknown fraction. A low-coverage trait result stays visible as partial evidence.

Use cultured oxygen phenotypes from BacDive, then explicitly labelled genomic predictions from metaTraits/porTraits and compatible trait resources. Strain-level evidence overrides a species-level generalization only when the strain is actually resolved. Deduplicate observations reused by multiple databases. The new card is labelled **Aerotolerance balance** inside the existing gut-community section, with an explanation of its relation to gut redox research. It is not a direct oxygen concentration or permeability test.

### 6.4 Functional redundancy and resilience

Add a functional redundancy table: for each validated route, count distinct resolved carriers and report their abundance-weighted effective number, coverage, and community-only unassigned evidence. For carrier contributions `w` normalized within a function, use `exp(-sum(w_i ln w_i))`. This is redundancy among resolved contributors, not a guarantee of ecological resilience.

Add an independently versioned research **resilience-associated community profile** using the published HACK resources where source-compatible inputs are available [E10]. For a sample-level output, implement the paper's **HACK-top-17 sample score**, not its different taxon-level HACK index. Pin the exact 17-taxon manifest and native rank calculation from the released code; reproduce its examples, then freeze the external cohort used for sample abundance ranks. Never rank against the current run's samples. Show the resulting research score and an eligible reference percentile if available, with coverage. It remains separate from existing scores. Do not manufacture a single-sample recovery probability by averaging diversity and favourable-looking species.

Actual recovery is a longitudinal feature. For a recorded perturbation with at least three eligible baseline samples, define baseline `b` as the arithmetic mean composition followed by closure. Define the observed baseline envelope `B = max_i BrayCurtis(baseline_i, b)`, enlarged to a separately validated technical-noise bound if one exists. This is an empirical envelope, not a 95% confidence interval from three samples. Compute `D(t) = BrayCurtis(p_t, b)`. Let `D_peak` be the maximum observed post-event distance and `t_peak` the earliest time attaining it. Only if `D_peak > B` and `D_peak > 0`, and only for `t >= t_peak`, calculate:

```text
observed_return_percent(t) = 100 * clip((D_peak - D(t)) / D_peak, 0, 1)
```

Label this **return toward the pre-event community**. Earlier post-event points show displacement only. This is a retrospective trajectory: additional samples can revise the observed peak and the new trajectory, while previously issued reports remain immutable. It does not establish that the old baseline was healthy. With fewer samples, still show changes and similarities; do not suppress the longitudinal section simply because recovery cannot yet be estimated.

### 6.5 Longitudinal views — A13

Support the last three measurements on compact metric cards and all available timepoints in an expanded chart. Keep original historical results immutable. Record comparability as `native`, `uniformly_reprocessed`, `bridge_validated`, or `descriptive_only`.

The fingerprint includes specimen/participant ID, collection date, stool form, extraction/library methods, storage, host filtering, read depth, tool/database/model versions, and recorded diet/antibiotic/supplement events. Repeated processing of identical FASTQs must retain the same specimen identity and must not be treated as additional biological timepoints.

For comparable abundance vectors:

```text
BrayCurtis(p,q) = sum(abs(p_i-q_i)) / sum(p_i+q_i)
Jaccard_distance = 1 - size(P intersect Q) / size(P union Q)
```

Use a fixed, documented detection rule for presence sets. Empty unions are undefined, not perfect agreement. Show raw changes separately from reference-percentile changes. A new model/reference version starts a labelled comparison segment unless samples have been uniformly reprocessed in an independent analysis lane. Reprocessing must not overwrite prior reports.

Follow-up suggestions explain what the user wants to learn and what exposure changed. Implement separate templates for routine review, a new probiotic/prebiotic, a major dietary change, and recent antibiotics. Do not hardcode one universal retesting interval as though it were established for every organism, intervention or symptom.

### 6.6 Larger reference and interpretation resources

Use curatedMetagenomicData 3 [E08] to discover eligible cohorts and metadata for **new** reference measurements. Its released profiles are not interchangeable with every profiler version; use compatible native outputs or reprocess public raw reads for the new lane. Preserve participant/study grouping and age/body-site eligibility.

Add the 2025/2026 large-scale health/diet association resource [E09] to organism explanations and action context. Its published 661-SGB ranking and public profiles can inform a separately labelled association view after taxonomy mapping. An association ranking does not make each organism intrinsically good/bad, and it does not replace any existing score. Resolve SGB splits and missing mappings explicitly rather than moving ranks between similarly named organisms.

### 6.7 Host-DNA context

Add a new host-DNA provenance card from trusted pre-filter QC counts or an eligible raw-read branch [E11]:

```text
host_DNA_percent = 100 * qualifying_human_aligned_reads / all_qualifying_input_reads
```

Use the same counting unit in numerator and denominator. If only host-filtered reads are available, state `upstream_host_filtered` unless original counts can be imported. Do not infer original host fraction from residual reads. Only aggregate counts enter the report; existing host-sequence handling remains unchanged. Explain the research association with inflammation without treating the fraction as a substitute for calprotectin or inferring cellular origin from ordinary sequencing.

## 7. Organism, strain, assay-coverage and resistance additions

### 7.1 Complete organism explorer — A10

Keep the existing organism pages and their values. Add a searchable, sortable view that joins existing lanes through accession/cluster mappings, not string similarity alone. Each row contains accepted name; original label; synonyms; stable taxonomy IDs and release; full lineage; rank; source lane; raw abundance and denominator; confidence/coverage; strain/subtype evidence; interpretation links; and detail-page link.

Genus and family rollups must disclose their membership and exclusions. A genus total cannot be combined with its own species rows as independent components of one pie chart. A subspecies such as *B. longum* subsp. *infantis* is nested within its parent when the method supports that hierarchy. GTDB suffixes and SGB placeholders are preserved as identifiers; resolving a display name is not a new detection. Species splits/merges require an explicit one-to-many mapping state.

Show all supported organisms, including low-abundance and unresolved ones, with filters rather than a silent row limit. Add explanatory cards for organisms lacking prose using curated traits and cited associations. `Unknown` is a meaningful interpretation category, not evidence of harm. Clinical effects and probiotic effects may be strain-specific. A species label alone does not establish a commercial probiotic strain, toxin carriage, morphology, or a disease-causing subtype.

Add a role-composition/interpretation-coverage bar over one declared abundance vector, using the unchanged existing classifications with their existing meanings. Show favourable-associated, conditional/context-dependent, opportunist/concern and unknown categories only where the existing registry supports the mapping; retain unassigned mass. Do not reclassify organisms to make a larger green segment.

Also provide a closed **34-name sentinel browse view**: *Akkermansia muciniphila*; *Bifidobacterium*; *B. adolescentis*; *B. longum*; *Blautia*; *B. obeum*; *B. wexlerae*; *Faecalibacterium*; *F. prausnitzii*; *Lactobacillus*; *Roseburia*; *Ruminococcus*; *R. bromii*; *Agathobacter rectalis*/*Eubacterium rectale*; *Alistipes finegoldii*; *A. putredinis*; *Anaerobutyricum hallii*/*Eubacterium hallii*; *Bacteroides*; *B. fragilis*; *B. thetaiotaomicron*; *Oxalobacter formigenes*; *Phocaeicola vulgatus*; *Prevotella*; *Prevotella/Segatella copri*; *Streptococcus thermophilus*; *Enterococcus faecalis*; *Escherichia coli*; *Klebsiella pneumoniae*; *Ruminococcus/Mediterraneibacter gnavus*; *Bilophila wadsworthia*; *Desulfovibrio*; *Fusobacterium nucleatum*; Proteobacteria/Pseudomonadota; and *Methanobrevibacter smithii*. These require searchable identity/status views, not a positive detection in every sample. Expand abbreviated genus names to full canonical names in the registry.

Add a **taxonomy-resolution view** that accepts organism names and resolves them to the current catalogue: exact identity, accepted synonym, cluster-level correspondence, split/merge ambiguity, supported detection, below quantification limit, no supported detection, or not assessable. Use every entry in the section-16.4 taxonomy fixture to test name resolution. It is not a list of organisms that must be found in every sample.

The current detection lanes remain operational. Validate each additional detection method before exposing its extra calls: use mocks, close relatives, dilution series, read-length/depth matching, and negative controls. A union of classifier names is not automatically a union of true species. Existing counts and scoring remain unchanged; an additional evidence-qualified inventory can expose independently supported extra detections without feeding them into old models.

Useful benchmark resources include the Sylph study/mocks [D01], CAMI [D02], and published shotgun gastroenteritis validation [D03]. Optional newer candidate profilers can be evaluated in a separate research lane; installing more tools is not itself an acceptance criterion.

### 7.2 Explicit named-target coverage

Build an extension view over the current 486-target registry and any verified additions. A target absent from PDF prose is not necessarily absent from the database. First resolve aliases and inspect installed references, assay compatibility and actual execution. Every requested name gets a row with `available_reference`, `analysis_ran`, `assay_compatible`, `evidence_state`, and `reason`.

At minimum, the following named targets/groups must be discoverable in the new view, without silently treating every listed organism as pathogenic:

| Group | Required searchable names and distinctions |
|---|---|
| Enteric/opportunistic bacteria | Enterobacteriaceae; *Klebsiella*, *K. pneumoniae*, *K. oxytoca*; *Escherichia coli*; Shigella/*E. flexneri* and Shigella/*E. dysenteriae* taxonomic aliases; *Salmonella enterica* and serovar Typhi; *Citrobacter*, *Enterobacter*, *Morganella*, *Raoultella*; *Campylobacter* and *C. jejuni*; *Helicobacter pylori*; *Vibrio* and *V. cholerae*; *Yersinia enterocolitica*, *Y. pseudotuberculosis*, *Y. pestis* |
| Additional opportunistic/context bacteria | *Streptococcus*, *Staphylococcus*, *S. aureus*; *Pseudomonas aeruginosa*; *Haemophilus influenzae*, *H. parainfluenzae*; *Enterococcus faecium*, *E. faecalis*; *Clostridioides difficile*; *Clostridium perfringens*; *Acinetobacter baumannii* |
| Fungi | *Candida* and *C. albicans*; *Aspergillus*, *Cryptococcus*, *Saccharomyces*, *Rhodotorula*, *Saprochaete*, *Malassezia*, *Geotrichum*, *Microsporum*, *Trichophyton*; Microsporidia with recognized human species where the assay covers them |
| Protozoa | *Blastocystis*; *Cryptosporidium*; *Entamoeba histolytica*, *E. dispar*, *E. coli*, *E. hartmanni*; *Giardia*; *Cyclospora cayetanensis*; *Balantidium/Balantioides coli*; *Chilomastix mesnili*; *Dientamoeba fragilis*; *Pentatrichomonas hominis* |
| Helminths | *Strongyloides stercoralis*; *Taenia* and *T. solium*; *Schistosoma*; *Fasciola*; *Hymenolepis*; *Dipylidium caninum*; *Diphyllobothrium/Dibothriocephalus latus*; *Enterobius vermicularis*; *Mansonella*; *Ancylostoma duodenale*; *Ascaris lumbricoides*; *Necator americanus*; *Trichuris trichiura* |
| DNA viruses | Enteric adenovirus F40/F41; cytomegalovirus; Epstein–Barr virus, each with specimen-context interpretation |
| RNA viruses | Rotavirus A; astrovirus; norovirus GI/GII; sapovirus — require an RNA/cDNA-capable input or imported assay for an assessed result |

Normalize misspelled input aliases, including `Yersinia enterolytica`, to the correct accepted identity while preserving the input string in audit metadata. “Microsporidium spp.” is not sufficient species resolution for all Microsporidia. Some entries are skin-associated, nonpathogenic colonizers, or poorly suited to stool detection; preserve that context. A DNA-only sample must not acquire a negative RNA-virus result through a missing-read default.

Support genus-level `Cyclospora` lookup as well as *C. cayetanensis*, and the historical name *Diphyllobothrium latum* alongside *Dibothriocephalus latus*.

The extension's assay table is a coverage and evidence view. It neither clears an FMT donor nor replaces the current donor-matching algorithm. Retain clinical context for low-depth organism detections and distinguish presence, subtype evidence and host disease.

### 7.3 Resistance summary — A11

Read existing AMR calls, then add a transparent class/mechanism summary using source ontology identifiers [D04–D06]. Include beta-lactamase classes A/B/C/D, macrolide–lincosamide–streptogramin B, fluoroquinolones, ciprofloxacin as a subordinate drug view, vancomycin, sulfonamide/trimethoprim, polymyxin/colistin, and aminoglycosides. Enzyme class, antibiotic class and individual drug are different levels of the hierarchy.

For every class show:

- Unique qualifying determinant count and sequence-family identity.
- Fragment-supported abundance in one declared unit and denominator.
- Acquired gene versus validated resistance-associated mutation versus unresolved homolog.
- Mechanism and drug specificity from the curated source.
- Carrier/linkage evidence and unresolved fraction.
- Assessed/not-assessed status and database/tool release.

Deduplicate alleles and parent/child ontology terms for aggregate richness. Do not count ciprofloxacin again as an independent class when it is already in the fluoroquinolone total. Use `AMRFinderPlus` with its supported assembled nucleotide/protein inputs; it is not a direct FASTQ classifier. Where using read mapping, use an appropriate read-based workflow and its own validation. Wild-type housekeeping hits or species identity are not sufficient evidence of phenotypic antibiotic resistance.

An optional absolute-load import can add copies per gram when a suitable laboratory assay or spike-in supports it. Do not convert relative sequence abundance into absolute organism/AMR load without that information.

### 7.4 Previously unfilled analysis slots

Inspect the repository's existing KO, CAZyme, strain and reference-dependent slots and implement the missing annotation capabilities required by this specification. Add a versioned general functional table and new evidence pages from compatible HUMAnN/orthology analyses where needed. Preserve the old scores and result objects: new observations appear as companion evidence unless a future, separately authorized model update explicitly consumes them.

For unresolved strains, show the strongest supported level: species, clade, dominant-population consensus, genomic linkage or named reference strain. Add targeted reference/assembly analyses where they can resolve useful new evidence. Do not turn insufficient depth into a named strain call, and do not remove already supported strain outputs.

## 8. Consolidated actions, foods and products — A12

### 8.1 Preserve and extend the evidence library

Retain the current intervention system and its options. Add every existing and new supported finding as a typed trigger, including organism, pathogen, biofilm, mycobiome, function and symptom/context findings, and add a consolidated planner over their cards. Section 19 defines the binding all-evidence retrieval and prioritization contract. **Human trials are not a prerequisite for a research option to appear.** Include human, animal, ex-vivo, in-vitro and mechanistic evidence, each accurately identified. Absence of a clinical trial must not become “no possible intervention.”

An option remains specific to the actual compound, food preparation, strain, substrate, formulation, model and endpoint studied. An in-vitro inhibitory concentration is not an oral dose; a chemically modified starch is not every resistant starch; a combined probiotic result is not automatically an effect of either strain alone. Preserve conflicting and null evidence beside positive evidence.

Each new evidence edge must contain:

```json
{
  "edge_id": "stable-source-target-intervention-id",
  "action_entity_id": "exact-ingredient-or-strain-id",
  "target_id": "canonical-taxon-reaction-symptom-or-lab-id",
  "target_direction": "increase|decrease|change|null|mixed",
  "endpoint": "what was actually measured",
  "study_model": "human_trial|human_observation|animal|ex_vivo|organoid|in_vitro|biochemical|in_silico",
  "population_or_model": "explicit study population/model",
  "formulation": "exact preparation or combination",
  "study_exposure": "source protocol, if verified",
  "effect_estimate": null,
  "uncertainty": null,
  "source_url": "stable primary publication URL",
  "source_verification": "full_text|abstract_only",
  "gut_delivery_evidence": "human_site_measurement|human_live_recovery|animal_site_measurement|direct_site_administration|oral_outcome_observed|strain_dna_recovery_only|simulated_digestion|acid_bile_assay|formulation_rationale|not_applicable_nonviable|route_mismatch|unknown|contradicted_for_tested_preparation",
  "delivery_source_urls": [],
  "tradeoff_edges": [],
  "applicability_notes": [],
  "display_status": "research_option"
}
```

Enum alternatives in this illustrative record describe allowed values; use one value per field in production. Store secondary reviews separately as literature-navigation records, not primary experiment edges.

### 8.2 Planner behavior

1. Collect all eligible options for the user's chosen goals and existing/new findings, with links to existing cards and per-finding search-coverage records. Deduplicate ingredients, strains, combinations and citations while retaining different formulations as distinct entities.
2. Present a concise **start here** list of one to three practical options, then a complete expandable catalogue. Sort deterministically using goal relevance, explicit contraindications, user preference/tolerance, evidence applicability, formulation match and feasibility. Expose the sorting reasons; this is prioritization, not a predicted probability of benefit.
3. Keep food, prebiotic, probiotic, supplement/compound, herb and lifestyle categories. Include clinician-review medication options already supported by the existing registry without changing their logic.
4. Show studied positive effects, possible undesired effects, and uncertainties together. Counts of favorable and unfavorable associations must not be interpreted as clinical net benefit.
5. Respect supplied allergies, medication interactions, age applicability and prior intolerance. Preserve the distinction between an unknown patient fact and an explicitly declared effective default under section 10. A reported inulin intolerance changes its ranking and suggests supported alternatives; it does not erase its measured genetic capacity or diagnose a disorder.
   Confirmed contraindicated/allergenic options remain visible as research evidence with the reason, but are not eligible for that person's new actionable “start here” list. This affects the new planner only and does not delete existing evidence records.
6. Include experimental antimicrobial candidates in the normal target-specific option list when supported by laboratory/animal evidence. Do not require a human trial or a verified retail formulation, or silently hide these candidates behind a human-only default. Describe actual targets and delivery evidence. Do not label a generic ingredient an eradication guarantee or remove an existing option solely for lacking an RCT.
7. Separate a study's exposure from a personal dosing instruction. A source protocol can be shown in an expandable research detail with its model and formulation. The planner must not automatically convert a cell-culture concentration or animal mg/kg exposure into a human regimen.
8. Attach a monitoring question to each proposed action: symptom tolerance, repeated compatible microbiome measurement, or a directly measured clinical endpoint. Do not promise colonization or a specific score change.

### 8.3 Food checklist and product guide

Generate a deduplicated shopping/food checklist from the chosen actions. Group by plant foods, whole grains, legumes, nuts/seeds, fermented foods, named fibre preparations and other relevant categories. Each item links to its target findings, study evidence and alternatives. Merge cocoa/cacao aliases and legumes/beans overlaps without merging materially different products such as cocoa husk feed and cocoa powder.

A probiotic/product record includes exact strain designation, taxonomic synonyms, tested combination, formulation, viable versus inactivated preparation, available quality evidence, and commercial link separated from the research citation. A retail species name without a strain match cannot inherit a strain-specific clinical claim. Every efficacy claim must link to evidence for the exact strain, formulation and endpoint. Product matching is optional; ingredient/strain-level guidance must work without commercial links.

The intervention vocabulary appendix defines the required food and ingredient concepts for the action planner. The original **60-record evidence seed matrix**, **17 SYN records**, and **85 indexed target-specific TP/BP/TM/BU records** are embedded later. Repeated papers or experiments across these registers must deduplicate at observation level, while preserving all index aliases. Its records were verified against publication records and complete abstracts; entries marked abstract-level require full-text extraction before adding exact strain effects, concentrations, effect sizes or directions not supported by the abstract. Supported abstract-level facts can populate explicitly labelled research cards immediately.

### 8.4 Additional current intervention research

Add these primary sources to the planner and its model-development fixtures:

| Source | What it adds | Implementation use |
|---|---|---|
| Deehan 2020 [R01] | Different resistant-starch structures produce different microbial/SCFA responses | Preserve substrate/formulation specificity |
| Wastyk 2021 [R02] | Human fibre and fermented-food interventions with microbiome/immune endpoints | Evidence-backed dietary-pattern cards |
| Lancaster 2022 [R03] | Individual and dose-dependent effects of distinct isolated fibres | Model multiple responses rather than all-fibre equivalence |
| Sardá 2025 [R04] | Baseline-community-dependent responses to unripe banana flour and inulin | Show responder uncertainty; do not copy PICRUSt predictions as measured shotgun functions |
| Song 2025 [R05] | Fibre-response prediction in prediabetes with sequence accessions | Reproduction research with phenotype/access restrictions retained |
| Connors 2026 [R06] | Experimental fibre–microbe combinations and model-guided community functions | Optional scenario research, not an automatic human treatment prediction |

These supplement, rather than replace, the existing treatment-evidence library.

### 8.5 Personalized synbiotic options

Implement the section-11 candidate engine as part of this planner. Include exact strain–substrate combinations, relevant component-only alternatives and transparent ranking tied to the person's findings. Use the embedded SYN evidence seeds and S-source assets. The core evidence mode and the broad target-specific expansion in section 19 are required; metabolic simulation is an optional supporting layer. Integrate its output into report section 15 and the existing metabolism details, with no additional top-level report section.

## 9. Host-signaling and symptom/condition navigation

Add a context explorer that links to existing scores and the new measured components. Each card separates sample observations, literature mechanisms, missing information and possible follow-up. A marker count is a navigation count, not a new disease probability.

### 9.1 GLP-1-related microbial mechanisms — A07

Include separate channels for SCFA-related signaling, bile-acid transformations, indole-related effects, and experimentally supported microbial proteins. Link the existing relevant scores unchanged. Host FFAR2/FFAR3 and other endocrine receptors are not microbial FASTQ targets. Indole effects depend on exposure and context [F43]; do not assign every indole route a universal GLP-1-stimulating direction.

Use Tolhurst's receptor experiments, Chambers' targeted colonic propionate intervention [C01], and the P9 experimental protein study [C02] as distinct mechanisms. A specific protein detector requires the experimentally tested sequence and near-homolog validation; a species abundance is not a substitute. Report supported gene/protein **capacity**, not expression, secretion, host GLP-1 concentration or equivalence to GLP-1 medication.

The panel must display quantitative component readings when available, plus evidence coverage. Do not collapse unlike, overlapping pathways into an invented hormone concentration. Optional measured GLP-1 imports must specify active/total assay, specimen, fasting/postprandial timing, collection handling and units.

**Required experimental P9 sequence panel:** implement `host_signaling.p9_sequence_potential` with the following verified anchors [C02]:

| Field | Exact value |
|---|---|
| Native locus | `Amuc_1631` |
| Protein | UniProt `B2UM07`; GenBank `ACD05451.1`; 748 amino acids |
| Genome/CDS | `CP001071.1`, `complement(1965361..1967607)`; 2,247 nt including stop |
| Original primary supplement | Yoon 2021, `Supp Table 3a`, Excel row 26; expression primers in `Supp Table 6`, row 11 |
| Protein sequence URL | `https://rest.uniprot.org/uniprotkb/B2UM07.fasta` |
| Genome record URL | `https://www.ebi.ac.uk/ena/browser/api/embl/CP001071.1?download=true` |
| Protein SHA-256 | `8618a2f5066940d6ae5bd5f99859f478a6dbc3036ec3b77eb59cc3412160f6a2` |
| CDS SHA-256 | `fb313b8d177bed36da19ae85884357114f00e29cf47387018b3d39fc55f95bdc` |
| Hash normalization | Uppercase unwrapped sequence, no newline; CDS in its coding orientation |
| Wrong-identifier negative control | `Amuc_1831`, `ACD05649.1` / `B2UN38`, 98 aa — a different protein, not an accepted P9 alias |

Quantify gene support, breadth, diagnostic sequence evidence and competitive mapping against related S41 proteases. Divergent homologs remain unresolved unless function-specific evidence supports transfer of the annotation. Report experimental P9 genetic potential; neither protein secretion nor the person's GLP-1 response is measured by these reads.

### 9.2 Endogenous ethanol and auto-brewery context

Add a bacterial/fungal fermentation panel using the 2026 observational cohort [C03], the high-alcohol-producing *K. pneumoniae* experiments [C04], and human microbial-ethanol/liver evidence [C05]. Include mixed-acid fermentation, heterolactic fermentation, ethanolamine-associated pathways, ethanol formation/consumption and acetaldehyde-related transformations where specific routes can be resolved.

Public research assets include **PRJNA1257465** and `https://github.com/hlfreund/ABS_CHsu`. The repository contains `adhE` mapping scripts but references an external `bacteria_adhE_gene_labeled_seqs.fna`; resolve its actual sequence accessions before claiming the paper-specific reference panel is installed. Public code is not proof every input asset is bundled.

Generic `adhE` or fungal presence supports neither a high-alcohol strain identity nor a diagnosis of auto-brewery syndrome. Display route evidence, resolved carrier evidence, and relevant optional measured ethanol/clinical context. Preserve existing liver and fungal scores.

### 9.3 Complete context inventory

| Context | Existing/new observations to link | Required distinction and primary support |
|---|---|---|
| Anxiety and depression | Existing depression, GABA, indole/tryptophan; new pathway detail | Human associations versus experimental behavior; no stool-derived anxiety diagnosis [F15, C06] |
| Autism research | Existing relevant observations; B1/quinone pathways and multikingdom data | Pediatric study population; diet and metadata access; no unverified adult classifier [C07–C08] |
| Serotonin, melatonin and sleep | Direct microbial synthesis evidence, host modulation, tryptamine pathways, optional symptom diary | These are distinct molecules/mechanisms, not brain concentrations [F17–F18, C09–C11] |
| Atherosclerosis | Existing cardiovascular profile, TMA and related functions | Link existing disease score unchanged |
| Hypertension | Existing profile and observed functional context | Reuse existing result, not a duplicate score |
| Auto-brewery/endogenous ethanol | New ethanol panel and existing fungi/strain evidence | Genetic potential versus measured ethanol/clinical state [C03–C05] |
| Constipation | Existing IBS-C, methane and new substrate/utilization data | Stool function versus breath testing and clinical symptoms |
| Diarrhea | Existing IBS-D/pathogen findings, bile transformations and stool form | Bile-acid-associated subset versus all diarrhea [C12] |
| Abdominal pain | Existing histamine/IBS, new fermentation context | Strain-dependent experimental mechanisms, not generic species causation [C13] |
| Wheat/fructan sensitivity | Existing celiac/IBS plus fructan utilization, symptom history | Celiac disease, wheat allergy and non-celiac sensitivity remain separate [C14] |
| IBS | All existing IBS/IBS-C/IBS-D/IBS-M scores | These are already present; add navigation and substrate context only |
| Thyroid context | Published overlapping observations; optional TSH/free T4/antibodies | Host blood/transcript features cannot be inferred from stool DNA; separate thyroid conditions [C15–C16] |
| Food allergy/immune tolerance | Existing eczema/urticaria and SCFA results, relevant organisms | Mouse transfer/consortium evidence does not establish personal allergen tolerance [C17–C19] |
| GLP-1-related physiology | Section 9.1 | Host endocrine response is distinct from microbial genetic potential |
| Eczema | Existing skin profile and observed functions | Preserve existing score; add linked evidence/actions |

The 2024 autism resource includes **PRJNA943687** and `https://github.com/qsu123/ASD_multi-kingdom_diagnosis`. Preserve the 2024 metadata-access correction and August 2026 disclosure correction in provenance. Public reads alone are not a frozen deployable classifier. This release's default deliverable is an evidence/context card; reproducing a new classifier is a separately labelled research capability, not an alteration of existing scores.

## 10. Your information, inputs and assumptions — A14

### 10.1 One combined report section, with useful downstream connections

Implement **report section 25: Your Information & Report Context** near the end of the document. Move the existing **“What this report knew about you, and what it assumed”** content into this section. Extend that same input register to support optional questionnaire answers and laboratory measurements. There is no separate top-level laboratory-results section.

The previous report section 25, **“What this report cannot tell you,”** becomes section 26. Subsequent existing sections shift by one: **“Sequencing quality and reference coverage”** becomes 27; **“How this test performs”** becomes 28; **“Technical data”** becomes 29. These are report section numbers, not the numbered chapters of this specification. Use stable semantic anchors and generate display numbering from the final document structure.

The combined section begins with a compact summary of supplied information, effective assumptions and inputs not supplied, followed by a readable table: **Information · Value used · Supplied or assumed · What it affects**. Every influential assumption must be visible. Group optional missing inputs into concise categories; do not pad the report with pages of empty laboratory charts. A supporting expandable table or appendix within this same section can list all recognized fields.

The current application works without external laboratory data. That remains the default path. Imports are scaffolding with an actual purpose: when a relevant measurement is supplied, show it beside the associated finding, explain its relationship to that finding, and use it in a documented interpretation or action rule where supported. Section 25 is the central record and provenance view; the relevant existing domain section is where its practical significance appears.

### 10.2 Supported patient and sample context

Extend the existing sample manifest and metadata intake rather than building a competing source of patient facts. Support the following categories when supplied:

| Input category | Useful fields | How it is used |
|---|---|---|
| Participant and specimen | Stable participant/sample ID; collection date/time; sample type; supplied demographics; country/region; height and weight with dates where useful | Link records correctly, establish comparability, select eligible references through existing supported interfaces, or contextualize a new finding |
| Symptoms and goals | Stool form/frequency; constipation/diarrhea; bloating, pain and onset; patient priorities; recent symptom changes | Organize context cards, link symptoms to relevant evidence and define monitoring questions |
| Diet and tolerance | Broad dietary pattern; specific fibre intake or intolerance; fermented foods; usual restrictions; recent changes; known allergies | Explain substrate findings and rank practical alternatives without inventing exact intake |
| Medications and supplements | Exact drug/product/strain where known; start/stop dates; relevant timing; antibiotics; acid suppressants; metformin; laxatives/bowel preparation; probiotics/prebiotics; antimicrobials; relevant other medicines | Apply existing context handling, distinguish exposure effects and annotate action applicability |
| Clinical history | Supplied confirmed diagnoses and relevant procedures; pregnancy/breastfeeding where applicable; immune status and other existing clinical-context flags; family history when a supported interpretation uses it | Add declared context and route relevant actions; a microbiome pattern is not a confirmed diagnosis |
| Collection and processing | Stool collection/storage conditions if known; extraction/library identifiers; provider host filtering; technical replicates; supplied QC metadata | Assess comparability and explain which inputs were available |
| Longitudinal exposures | Dated diet changes, medications, supplements, procedures and other relevant perturbations | Annotate trends and distinguish a sample change from a method change |
| Optional laboratory results | The supported categories in section 10.4 | Link measured endpoints to genetic potential and relevant existing findings |

Use progressive optional questions: ask only for fields relevant to the chosen analysis or action, permit skipping, and preserve the ability to generate a full report from the current inputs. Do not collect unrelated personal information merely because a flexible schema can store it.

### 10.3 Effective defaults and the known-versus-assumed contract

1. Resolve an explicit, applicable supplied fact first, using specimen date and participant identity. Conflicting supplied records require a documented resolution or a visible conflict; do not overwrite one silently.
2. Reuse the production application's existing defaults and metadata rules for existing capabilities. Import the existing assumptions register into section 25; this specification does not silently redesign its rules.
3. For a missing input that influences a **new** capability, define a reasonable versioned fallback in that capability's input contract. State the effective value, rationale and exactly which result or action uses it. A missing optional answer must not stop the whole report.
4. Default examples: a new longitudinal view uses the available dated samples; a new dietary scenario uses its declared reference medium when intake is absent; a new reference comparison uses the documented eligible reference population when tighter matching information is absent; a model may use its frozen training-set imputation rule when that rule was part of its validated workflow. Where several plausible assumptions materially change a new result, show the scenario range or a sensitivity note.
5. Keep the recorded patient fact separate from the computational fallback. “Assumed for this calculation” must never become “patient confirmed.” Do not invent a precise demographic value, diagnosis, allergy history, medication history or intake amount simply to fill a field. Existing defaults remain visibly identified as defaults.
6. Laboratory measurements have **no imputed normal value**. If absent, display **“Not supplied”** in section 25, retain null in the measurement record, and continue the sequencing analysis. Do not generate normal/abnormal classifications for absent assays.
7. When a fact replaces an assumption, use the application's existing metadata-consumption behavior and declared dependencies. Recompute only affected outputs through supported interfaces; record the change and retain previous reports. New laboratory imports do not become undocumented inputs to an existing score.

Each context field uses this logical record, adapted to the existing schema:

```json
{
  "input_id": "diet.scenario_medium",
  "category": "diet_context",
  "label": "Dietary medium used for the optional simulation",
  "participant_id": "SAMPLE_OWNER",
  "sample_id": "SAMPLE",
  "reported_value": null,
  "effective_value": "declared_reference_medium",
  "unit": null,
  "input_state": "assumed_default",
  "source": "scenario_configuration",
  "recorded_at": null,
  "valid_at": null,
  "default_id": "scenario.default_medium",
  "default_version": "resolved_from_installed_model",
  "default_rationale": "No intake data supplied; use the model's declared reference scenario.",
  "derivation": null,
  "affects": ["scenario.substrate_response"],
  "linked_results": [],
  "limitations": ["This describes a modelling assumption, not the person's actual diet."]
}
```

`input_state` supports `supplied`, `derived`, `inherited_default`, `assumed_default`, `not_supplied`, `not_applicable` and `conflicting`. Supplied values, derivations, defaults and conflicts have independent provenance. Model/configuration versions in the example must resolve to actual installed identifiers in production.

Every influence declaration contains its target metric/action, the consuming rule, whether it changes a calculation or interpretation or merely adds context, and a reader-facing explanation. If a stored field has no implemented consumer, label it **“Recorded; not used in this report”**. Do not imply integration merely because a value was imported.

### 10.4 Optional measured inputs and their destinations

| Input panel | Supported imported analytes | Integration in the existing report |
|---|---|---|
| SCFAs | Acetate, propionate, butyrate, valerate and total SCFA concentrations; measured fractions when supplied | Fermentation cards in report sections 10/20, with DNA capacity and measured concentration as separately labelled values |
| Inflammation/immune markers | Fecal calprotectin, lactoferrin, lysozyme, secretory IgA | Relevant existing condition/barrier context and associated detail cards; never an automatic disease diagnosis |
| Digestion | Fecal fat, pancreatic elastase, reducing carbohydrates/carbohydrate assay | Relevant digestion/substrate or symptom context, with direct-assay interpretation distinguished from microbial capability |
| Other stool chemistry | Occult blood/FIT or guaiac method, stool pH, measured beta-glucuronidase activity | The applicable condition, fermentation or enzyme-activity detail; preserve the assay-specific meaning |
| Targeted metabolomics | Bile-acid species, urolithin/equol challenge results, relevant aromatic metabolites, measured ethanol/acetaldehyde where available | Corresponding metabolic cards and their action/context explanations |
| Endocrine/nutrient context | GLP-1 with assay/timing, thyroid tests, clinically obtained nutrient markers | Corresponding signaling, vitamin or condition-context card; a microbial gene is not a hormone or nutrient concentration |
| Absolute microbial load | Validated qPCR, flow-cytometry or spike-in-supported load with denominator and units | Organism/ecology or AMR detail where units, timing and method permit a supported comparison |
| Molecular pathogen assays | Qualitative/quantitative PCR or other supplied DNA/RNA assay, with named target and specimen | Existing pathogen finding and coverage tables, clearly attributed to the external assay |

Every laboratory record retains analyte ID, original name, value and censoring (`<`/`>`) or qualitative result (`detected`, `not_detected`, `indeterminate`), original/normalized units, method, specimen, collection date, laboratory interval, source identity and document hash. Molecular records also retain target, nucleic-acid type, optional Ct and detection limit. Check participant/specimen identity, date compatibility and matrix before linking. Preserve original values; unit conversion requires analyte- and matrix-specific rules. Wet/dry stool weight, relative percentages and concentrations are not interchangeable.

For each supported input, implement a binding to the appropriate canonical finding, a comparison rule or a clearly stated context-only role, and bidirectional links between the finding and its section-25 record. Where an evidence-backed action rule uses the measurement, expose that rule and source. A supplied result may support, differ from, or be incomparable with genetic potential; explain rather than force agreement. Incompatible or historical measurements remain labelled as such.

In an ordinary report with no external tests, show **“External laboratory results: Not supplied”** and a compact category list in section 25. Avoid repeated empty charts or a generic sales pitch for more testing throughout the report. A relevant detail card may explain which optional measurement would clarify its endpoint. Importing data is optional and must never be required to render existing results.

Vaginal microbiome community state and other body-site measurements are not derived from stool. Future imports retain their real specimen and cannot be silently treated as gut measurements.

## 11. Personalized synbiotic options and model scenarios — A16

### 11.1 Required capability and placement

Build a **personalized synbiotic candidate recommender** from the person's existing findings, measured strain/function evidence, goals, supplied context and the study-specific intervention registry. The core evidence-matching mode is required and works without a metabolic solver. An optional model-assisted mode adds reproducible counterfactual comparisons. This is an experimental recommendation system, not an individually validated prediction of benefit.

Place **Personalized Synbiotic Options** inside existing report section 15, **What you can do**. Put substrate use, cross-feeding, model outputs and source details within report section 20, **Microbial Functions & Metabolism**, with reciprocal links. Do not create a new top-level report section. Integrate candidates with the existing action list and food/product checklist instead of repeating recommendations.

Personalization means explaining why an exact strain, substrate or combination is relevant to this person's supported findings. It does not mean that detecting a low species identifies a deficiency, that a commercial strain is identical to the resident strain, or that adding a probiotic guarantees lasting colonization. Include relevant preclinical options with their actual evidence setting; lack of a clinical trial alone is not an exclusion criterion.

Use ISAPP terminology precisely [S01]. Store `design_intent = complementary | synergistic | unspecified` separately from `definition_status = candidate | host_benefit_supported | definition_criteria_met | unresolved`. Genetic compatibility or increased culture growth alone does not establish a health-benefiting synbiotic. A statistical interaction test is useful evidence but is not, by itself, the definition of a synergistic synbiotic. Use **candidate combination** where the requirements have not been established.

### 11.2 Canonical entities and evidence contract

Extend the application's existing intervention, organism and source registries. These are logical entities, not a prescribed directory structure:

| Entity | Required fields and semantics |
|---|---|
| Strain | Stable local ID; current and publication taxonomy; exact strain designation; verified collection aliases; BioSample and versioned assembly/sequence accessions when available; genome identity evidence; viable/inactivated state; source references. A species proxy has a separate ID and never acquires a branded-strain identity. |
| Substrate preparation | Chemical identity/linkages, degree-of-polymerization range when reported, purity, residual sugars, source/manufacturing treatment, mixture proportions, dose/concentration with units and setting. Whole foods, isolated polymers and model exchange proxies remain separate. Unreported fields are null. |
| Formulation | Exact constituent IDs and proportions, delivery/viability information, studied preparation, product/lot/date if verified; atomic study-bundle identity. A combination's result cannot be assigned to each component separately. |
| Observation | Study/arm, biological donor/participant, replicate type, organism or community, medium/pH/oxygen/time, intervention and comparator IDs, endpoint/unit, actual effect/statistics or qualitative finding, positive/null/negative/conflicting state, source table/page/figure, access and extraction status. |
| Evidence edge | Subject strain/community + preparation + endpoint; `measured_growth`, `genomic_prediction`, `culture_metabolism`, `animal_host`, `human_ecology`, or `human_clinical`; exact versus proxy mapping; tested context; contrary observations. Multiple papers from one cohort do not become independent replication. |
| Candidate | Components/formulation, findings and goals it addresses, known interactions and contraindications, studied versus novel combination, population applicability, model coverage, ranking explanation, evidence/source IDs and optional exact-product match. |
| Scenario | Complete sample/model/medium/solver fingerprints, perturbation, comparator, feasible status, objective and units, production/growth outputs, sensitivity configuration and coverage. |

For growth datasets distinguish `tested_no_growth`, `not_tested`, `genomic_pathway_absent`, `no_added_combination_benefit` and `contraindicated`. None is an alias of another. A null statistical result is not proof of zero effect. Preserve biological donors separately from technical replicate counts.

The minimum candidate object is:

```json
{
  "candidate_id": "synbiotic.candidate.example",
  "kind": "candidate_combination",
  "components": [
    {"entity_id": "strain.example", "identity_level": "strain"},
    {"entity_id": "substrate.example", "preparation_match": "exact"}
  ],
  "design_intent": "unspecified",
  "definition_status": "candidate",
  "combination_status": "studied_exact",
  "personalization_basis": ["supported_finding", "supplied_goal"],
  "finding_ids": [],
  "goal_ids": [],
  "evidence_ids": [],
  "source_ids": [],
  "applicability": "population_mismatch",
  "eligibility": "context_caution",
  "eligibility_reasons": [],
  "goal_coverage_percent": null,
  "pair_supported_goal_coverage_percent": null,
  "rank": null,
  "rank_scope": "evidence_matching",
  "rank_policy_version": "targeted-option-ranking/2",
  "scenario_ids": [],
  "model_state": "not_run",
  "known_tradeoffs": [],
  "unknowns": [],
  "monitoring_question": null
}
```

This is a schema example, not a recommendation for a supplied sample. Null means unavailable; it must not become zero response, zero risk or an invented certainty.

### 11.3 Candidate generation and personalization

1. **Resolve the inputs.** Use existing measured taxa, strain objects, new substrate/function capacities, relevant findings and context. Preserve namespaces and assay coverage. Low sequencing depth, a species-only detection or an unresolved genome match cannot become an exact-strain call. Record assumed diet separately from supplied diet.
2. **Define goals.** Use user-selected goals first, then supported existing findings with a versioned explicit goal-to-endpoint map. For example, a supported low-butyrate-capacity finding may retrieve butyrate-oriented options; measured low stool butyrate and low genetic capacity remain different triggers. Do not invent a new universal direction for metabolites, diversity or individual organisms. Conflicting goals remain visible.
3. **Retrieve candidates.** Include exact studied pairs/formulations, supported single ingredients, and explicitly novel combinations of supported components. A novel combination needs either a relevant complementary rationale or a strain–substrate utilization/cross-feeding edge; do not enumerate an unsupported Cartesian product. Every returned candidate must have a source-backed edge to a goal. Whole-trial bundles remain atomic unless separate component evidence exists.
4. **Attach identity and applicability.** Exact strain evidence, species-proxy evidence, population matches and preparation mismatches remain distinct. Infant/animal/ex-vivo evidence is searchable and visible; it cannot silently become an adult clinical outcome. Delivery evidence is separately graded as measured in humans, simulated digestion, in vitro acid/bile, formulation rationale or unknown.
5. **Include comparators.** Evaluate no addition/current baseline, substrate alone, probiotic alone and combination wherever data/model support them. A missing component-only study arm remains missing. Never manufacture factorial evidence. Fibre alone or no clear preferred candidate may be the answer.
6. **Check relevant tradeoffs.** Incorporate explicit allergies, prior intolerance, known interaction/contraindication evidence and preparation contents such as residual lactose. Modelled feeding of a concerning resident organism is a hypothesis requiring a compatible strain/pathway/model edge, not a blanket rule that all fibre promotes pathogens. A disease resemblance score alone cannot exclude a strain or imply that it causes that disease.
7. **Keep options usable.** An explicit contraindication prevents placement in the person's actionable start-here list, with the source and reason shown. Research-only organisms, unknown commercial availability and population mismatches remain visible as such. Unknowns must not be converted into contraindications or silently treated as absent.
8. **Attach follow-up.** Each option identifies a relevant tolerance/clinical/ecological endpoint and what additional information would reduce uncertainty. Lasting engraftment requires appropriately timed longitudinal exact-strain evidence; growth in a model and detection during ingestion are different observations.

### 11.4 Transparent ranking; no invented treatment-success percentage

Provide ordered options with an explanation and separate evidence labels. Use the existing action planner's person-specific preferences and exclusions; the following deterministic policy fills the new synbiotic behavior.

For goals \(j\) with nonnegative weights \(w_j\), default equal weights unless supplied, define \(c_j=1\) when the candidate has a source-supported favorable edge relevant to that goal at its disclosed evidence level, otherwise 0. Define \(p_j=1\) only when the exact combination itself, rather than either component alone, supports that endpoint. Report:

\[
C=100\frac{\sum_j w_jc_j}{\sum_jw_j},\qquad
P=100\frac{\sum_j w_jp_j}{\sum_jw_j}.
\]

For a single-component option, P is null (not applicable); C is still computed. If no eligible goals exist or the total weight is zero, both are null. Reject negative or nonfinite weights. These are **goal coverage** and **pair-supported goal coverage**, not probabilities, clinical benefit magnitudes, or health scores. A conflicted goal retains its favorable and contrary edges and is marked conflicted; the coverage display must not hide this distinction. For each covered goal show the evidence setting. Do not give a mouse result the appearance of a human clinical response by using the same percentage without that breakdown.

Apply the relevance-first `targeted-option-ranking/2` policy in section 19.5. All evidence settings are eligible: do not use a universal human-to-culture evidence ladder, or classify every nonhuman observation as ineligible for the normal finding-specific list. Exact target, actual endpoint, formulation, independent replication, delivery and explicit tradeoffs determine relevance. Keep the study setting visible. Use C and P as descriptive summaries alongside the goal-specific ranking vectors; neither enters the rank vector or means a success probability. The schema has no blanket study-setting ordinal.

Return the complete ordered candidate list and the one-to-three-item start-here view from section 8. An experimental option can be prioritized with its exact evidence label. Explicit contraindications affect start-here eligibility, not evidence retention; unknown exposure is displayed, not transformed into a contraindication. Mechanistic-only concepts remain identifiable. If no candidate has a relevant empirical effect, say so precisely and retain supported mechanistic candidates and the search-coverage state.

Keep `evidence_matching` and `model_assisted` ranks separate. Compare model scenarios numerically only when model/objective/units, comparator and goal definition match. Within that scope order by descending median goal-directed delta, descending minimum goal-directed delta, then stable ID; expose the full sensitivity grid. Conflicting goals need separate tables or a disclosed user-weighted objective. Model outputs never silently change evidence rank or manufacture observed effects. Unknown tradeoffs are not zero harm, and effect magnitudes from incompatible endpoints cannot be compared.

For a novel multi-component design, keep studied bundles atomic, deduplicate constituents, enumerate supported additions with a configurable search budget, and re-evaluate the whole combination. A union of individually promising pairs does not inherit their individual outcomes. Compare dropping each component and compare to the best single intervention. Search budget is an engineering setting, not a clinically established maximum number of strains. Unsupported combinations stay unranked research concepts.

### 11.5 Open strain–substrate data and seed ingestion

Use the embedded SYN seed records below and sources S01–S15. The minimum pipeline imports exact experimental identity, preparation, comparator and outcomes, then independently reconciles strain genomes through NCBI/BioSample/culture aliases. No commercial product match is required for the recommender to work.

**Primary scalable training data:** Arzamasov 2025 [S02] provides 3,083 genomes, 68 reconstructed utilization pathways and experimental validation across 30 strains and 38 phenotypes. The genome-wide matrix is predicted; the growth-table observations are measured. Import raw growth Table 16 and substrate/metadata/validation tables separately. Split model evaluation by strain/lineage and experiment; retain assay negatives. Do not count thousands of predicted rows as thousands of experimentally validated strains.

**Supporting sources:** dbCAN-PUL supplies experimentally grounded polysaccharide loci; BacDive supplies strain identity and measured/predicted traits; Probio-Ichnos supplies literature-indexed acid/bile/adhesion and related evidence, not a substrate-growth or clinical-response matrix. IPDB is an optional genome catalogue. Resolve each edge to its actual experimental strain and assay. Missing catalogue traits do not imply incapacity.

The seed table is a starting registry, not an exhaustive catalogue of all possible combinations. Records with incomplete preparation or identity extraction remain valid evidence cards but cannot claim exact formulation matching until resolved. Record the unresolved field and source needed; do not invent it or call the whole feature unavailable. Incorporate newly curated pairs through the same schema and tests.

### 11.6 Model-assisted comparisons and experimental designs

Implement three distinct adapters; never merge their training domains:

- **MICOM/AGORA counterfactuals [R07–R08, S10–S11]:** run supported species/community metabolic scenarios for this sample, reporting the fraction of eligible bacterial abundance represented. Pan-species/pan-genus models are labelled proxies. A model-specific strain prediction requires a compatible strain model; species-level recommendations can still be shown with their limitation.
- **Strain–substrate compatibility [S02–S07]:** rank experimentally supported or genomic-predicted utilization and cross-feeding edges. Presence of a single generic CAZyme is not proof of polymer breakdown/uptake; substrate structure and complete supported loci matter.
- **MiRNN culture-design research [R06, S12]:** reproduce within the published fixed strain/fibre system. The released model does not accept arbitrary patient metagenomes as a validated personalization input. Use it for mechanistic candidate evidence and culture-domain benchmarks unless a separately validated adapter is developed.

For MICOM, reject negative/nonfinite abundances and a zero-total input, merge true aliases once, and normalize baseline abundance mass to one. After constructing the perturbation, assert its total equals one within absolute tolerance 1e-9 **before** any post-perturbation normalization; only normalize floating-point residuals within that tolerance. An invalid recipe cannot be repaired silently by dividing its total away. Report mapping loss and threshold exclusions before renormalization. Do not mix bacterial model coverage with fungal, viral or all-read denominators. Unsupported taxa remain unmodelled; high abundance coverage does not prove complete functional coverage.

Baseline medium, host-absorption assumptions, model version, solver/tolerances and objective are required inputs. With unknown diet, use a declared finite panel of media assumptions rather than inventing an observed diet. Model substrate uptake and degradation according to actual model chemistry; do not give every organism free glucose and label that resistant-starch utilization.

For each supported pair and each scenario compute the four conditions \(F_{00},F_{10},F_{01},F_{11}\), where indices mean probiotic and substrate absent/present. Hold nonintervention constraints fixed. State whether adding substrate increases carbon supply or substitutes equal carbon; never compare those experiments as equivalent. For every compatible endpoint output:

\[
\Delta_{pair}=F_{11}-F_{00},\quad
\Delta_{over\ substrate}=F_{11}-F_{01},\quad
\Delta_{over\ probiotic}=F_{11}-F_{10},\quad
I=F_{11}-F_{10}-F_{01}+F_{00}.
\]

`I` is an **interaction contrast in the model/experiment**, not proof of clinical synergy. For a small software fixture, values 1, 3, 4, 8 give pair change 7, improvement over substrate 4, improvement over probiotic 5 and interaction 2. Compute no contrast if a required arm is missing. No-addition, component-only and negative-response results must remain visible.

Keep raw production flux and flux divided by community growth as separate metrics. Their units and rankings differ; divide only by finite positive growth. Predicted flux is not stool concentration, absorbed exposure, administered dose or symptom improvement. Sensitivity across predeclared media/parameters is shown as min/median/max and rank changes, not a calibrated confidence interval or probability. Nonoptimal/infeasible solver status is unavailable, not zero production. Report model coverage and reasons for any suppressed scenario individually; unrelated candidates still work.

### 11.7 Pinned public reproduction and numerical fixtures

Use R08's final 2026 publication and the pinned source assets in S10–S11. Recreate input/media artifacts and lock a compatible environment before claiming an end-to-end reproduction. The following is an explicitly corrected paper-protocol configuration, not blind execution of the inconsistent Study A notebook:

```yaml
micom_version: 0.37.0
model_asset: agora201_refseq216_species_1.qza
cutoff: 0.001
tradeoff: 0.99
strategy: none
normalize_input_abundance: true
resident_fraction: 0.8
added_fraction_per_species: 0.04
added_species:
  - Akkermansia muciniphila
  - Clostridium beijerinckii
  - Clostridium butyricum
  - Bifidobacterium longum
  - Anaerobutyricum hallii
merge_existing_and_added_species: sum
required_total_abundance: 1.0
primary_replay_objective: raw_butyrate_flux
secondary_replay_objective: butyrate_flux_per_community_growth
```

This consortium is a **species-model perturbation**, not an exact-strain product or personal dose. The published grid uses six substrates plus no substrate, with/without consortium. Frozen model exchange caps in mmol/gDW/h are: `EX_inulin_m=6.14`, `EX_pect_m=0.4`, `EX_strch1_m=16.65`, `EX_dextrin_m=30.28`, `EX_cellul_m=0.37`, `EX_arabinoxyl_m=4.09`. These are total caps, not additions to the old cap. The last two are cellulose/hemp and arabinoxylan/psyllium proxies, not complete foods. Preserve the original exchange metadata; a missing reaction is not a successful scenario.

Handle these audited source issues explicitly:

1. Study A's notebook uses resident `.95` plus `.04×5`, while the paper and Arivale recipe use `.8 + .04×5`. Use the normalized paper protocol above and identify the repair. Do not copy total abundance 1.15 into production.
2. Model `.qza` media are not bundled; import the pinned CSVs. Assert duplicate `reaction`/`reaction.1` fields agree and populate all required identifiers for new exchanges. The model library is a separate download.
3. Methods/prediction code uses a growth threshold `.01 h−1`, while one figure caption says `.001`; the Arivale notebook also references an undefined aggregate `MCMM_engraftment`. Preserve this discrepancy. Any threshold sensitivity is research output, not an engraftment probability.
4. Public input has 154 unique sample IDs whereas reported output has 156. Whole-cohort replication requires an explicit ID-intersection audit; do not fabricate the two inputs.
5. No complete environment/solver lock is supplied. Build one, record license/solver status, then run a real example. Updating MICOM requires a separate versioned validation. Published-output replay alone is not a solver reproduction.
6. S5 uses Strict OOXML; support its namespace or convert with documented provenance. A parser returning zero sheets must fail loudly. Its treatment matrices contain growth-normalized outputs, distinct from raw production sheets.

**Public-output replay fixture:** S5 `production_but`, sample string `22001612560041`, standard-diet candidates only (exclude `High Fiber`). Values below come from released predictions, not a new OpenBiota simulation:

| Candidate | Raw butyrate flux, mmol/gDW/h | Change from no addition |
|---|---:|---:|
| No substrate or consortium | 6.5509159065 | 0 |
| Psyllium proxy alone | 19.8677059714 | 13.3167900649 |
| Maltodextrin proxy alone | 19.6188886985 | 13.0679727920 |
| Starch proxy alone | 19.4780226726 | 12.9271067661 |
| Inulin plus consortium | 12.1235454393 | 5.5726295328 |

Assert these values within rounding tolerance `1e-8`. In the separate growth-normalized objective the same sample's baseline is `124.022749302492`; maltodextrin alone `378.92775471395` exceeds psyllium alone `377.04681995744096`. Assert both orders under their correctly named objectives. This is a ranker/parser fixture; do not display it as an outcome measured in the user or a clinical recommendation.

For independent method checks, R07 provides public experimental projects `PRJNA937304`, `PRJNA640404`, `PRJNA939256` and `PRJNA1033794`. Reproduce the appropriate experimental medium and model version. R07 used MICOM 0.32.5/AGORA1.03 and distinct tradeoff/media settings; it is not interchangeable with the R08 configuration. Some human cohorts/phenotypes remain request-only; the core recommender must not depend on them.

### 11.8 Validation, reporting and completion

Required validation separates **identity**, **evidence extraction**, **analytical/genomic compatibility**, **numerical model reproduction**, and **human-response prediction**. Passing one is not passing the others.

- Verify exact seed identities and comparator outcomes, including the negative GOS and XOS factorial examples. A null interaction must not disappear because a component had a favorable result.
- Import both observed and predicted growth datasets; hold out strains/lineages and studies when training a new predictor. Do not leak replicated wells, repeated participants or a strain under a synonym into held-out data.
- Reproduce a public-output fixture and a genuinely solved model fixture separately. For a newly fitted model, compare against simple substrate-only, baseline-composition and nonpersonalized alternatives on held-out appropriate-domain data. Report endpoint-specific error, calibration where applicable, coverage and failures.
- An optional prospective personalized-response evaluation needs independently recorded baseline, actual formulation/exposure and outcomes. The present data do not establish a clinical success probability for new personalized combinations. Do not fabricate such a field.
- Every visible option shows its exact components, why it was retrieved, observed versus predicted effects, evidence setting, whether the combination itself was tested, key tradeoffs, missing information and links to relevant function details. Explain whether fibre alone or the combination ranked higher and why.
- Mark a commercial option as an exact study match only after its identity/preparation check passes. Products with incomplete composition can have a clearly labelled product-information card linked to relevant ingredient-family research, but cannot inherit exact-strain, exact-substrate or whole-formulation efficacy. A brand match or label does not independently establish viable strain content or efficacy. Unavailable research strains remain labelled research candidates, with no automated manufacturing/administration instructions.
- The evidence mode is complete only when SYN records are normalized, candidates/ranks are computed from real inputs and the report integration works. Optional simulations have their own feature flags and exact readiness states. An unimplemented model is not a successful computation and must not block the evidence mode.

A separate optional food-DNA adapter can retain **MEDI** [E16] for recent dietary-exposure context. Use its final 2025 code/reference manifest. Food-DNA absence does not mean zero intake and detected reads do not measure grams eaten. This adapter stays in existing diet/function details.


## 12. Report layout and interaction — A15

### 12.1 Integrate additions into the existing sections

The report remains one coherent document. New implementation capabilities do not receive new top-level sections by default. Place each new reading with its related existing readings, using subheadings, additional cards and continuation pages where needed. Preserve calculations and scientific distinctions while improving display titles, grouping and navigation. A new top-level section is justified only for a genuinely distinct subject with no sensible existing home. For this release, the only required new top-level section is the combined input/context section 25 described in section 10 of this specification.

Rename the cover's **“Metabolite Production”** label to **“Microbial Functions & Metabolism.”** Apply that umbrella consistently to existing report section 10, **“Your Microbial Functions & Metabolism at a Glance,”** and section 20, **“Your Microbial Functions & Metabolism in Detail.”** The broader title covers synthesis, degradation, substrate use, transformations and microbial signaling mechanisms. Keep short explanatory copy: **“What your microbes can make, break down and transform.”** Stable section and metric IDs remain unchanged even when display labels change.

Within report sections 10/20, group existing and new readings together under coherent subheadings:

- **Fibre & Dietary Substrates:** the complete 14-substrate panel, relevant existing carbohydrate functions and associated food opportunities.
- **Fermentation & Cross-feeding:** existing SCFAs plus acetate, lactate, succinate, hydrogen production/consumption and related network views.
- **Vitamins & Nutrients:** the complete nine-vitamin overview and synthesis/salvage explanations.
- **Protein & Nitrogen Metabolism:** existing and new amino-acid, ammonia, aromatic-metabolite and BCFA functions.
- **Gut–Brain & Metabolic Signaling:** relevant existing neuroactive compounds, GABA breakdown and the component-based GLP-1 mechanism panel.
- **Plant-Compound Conversion:** polyphenols, urolithin, equol and glucosinolate conversion.
- **Mucus, Bile & Other Transformations:** relevant existing functions plus reaction-specific mucin, lipid-A, bile, urate, polyamine, glutathione, ethanol and drug-metabolism detail. Use smaller subsections as needed; avoid a single overfilled card.

These are reader-facing groupings, not permission to merge unlike biochemical quantities. Each canonical metric has one primary home; other contexts link to it rather than repeating full pages or counting it as another measurement.

| Capability | Existing at-a-glance home | Existing detail/action home and integration |
|---|---|---|
| A01 — substrates | Report 10 | Report 20 substrate group; food/fibre choices in report 15 and the relevant detail card |
| A02 — fermentation | Report 10 | Report 20 fermentation group; related organism groups link from reports 5/16 |
| A03 — vitamins | Report 10 | Report 20 vitamins group; existing and new vitamins share one dashboard |
| A04 — protein/nitrogen | Report 10 | Report 20 protein/nitrogen group |
| A05 — neuroactive metabolism | Report 10 | Report 20 signaling group; relevant condition cards in reports 11/21 link to canonical readings |
| A06 — plant-compound conversion | Report 10 | Report 20 plant-compound group; organism-guild context remains in reports 5/16 |
| A07 — GLP-1 mechanisms | Report 10 | Report 20 signaling group, with clear gene-potential versus host-response labels |
| A08 — mucin/lipid A/bile | Reports 10 and 5, according to whether the reading is a function or organism group | Reports 20/16; related biofilm context links from reports 12/22 without inventing another top-level subject |
| A09 — ecology | Report 4, **Your gut community** | Add subordinate community-type, diversity, dominance, ratios, aerotolerance and resilience cards here; trait/guild detail links to report 16 |
| A10 — organism explorer | Report 7 organism inventory; report 6 when a finding needs attention; report 8 for named-target coverage | Existing organism/strain/pathogen detail in reports 17/18/19, with targeted additions in reports 13/23 for fungi |
| A11 — resistance | Existing Pathogens section, report 8 | Report 19 pathogen/AMR detail; add class and carrier-evidence subheadings |
| A12 — actions | Report 15, **What you can do about it** | Existing per-reading action blocks; one consolidated food/action checklist and formulation comparison in this section |
| A13 — time and context | Small trend views within each relevant existing overview | Full metric history with its existing detail card; symptom/condition navigation inside reports 11/21 or 14/24 as appropriate, not a separate catch-all section |
| A14 — inputs and optional measurements | Relevant supplied measurements beside their associated existing readings | New report 25 holds the unified input/assumption register and provenance; no separate laboratory chapter |
| A15 — navigation | Page 3 major-section contents; existing cover and overview cards | Metric-to-detail and back links throughout; complete detailed index/search in the existing technical/data area |
| A16 — personalized synbiotics | Report 15 Personalized Synbiotic Options, integrated with existing actions | Report 20 substrate/cross-feeding mechanisms, evidence and optional model detail; reciprocal links, no new top-level chapter |
| Host-DNA and assay context | Compact sample-context summary in report 25 where useful | Existing sequencing-quality/reference section, renumbered to 27; retain its units and input limitations |

All existing major sections not named in this mapping remain in the document and in applicable navigation. Generate the full section list from the report's actual structure rather than using this mapping as a whitelist. Biological-condition panels, organism pages, fungal findings and biofilm findings keep their existing homes. Detailed gene tables belong in the existing technical section, now 29.

There is no separate “Expanded capabilities” page or second general overview. Add new values to the appropriate existing overview, allowing that section to span additional pages when necessary. Do not truncate results to fit a predetermined page count.

### 12.2 Page 3: concise reading guide and clickable contents

Page 3 remains **“How to read this report.”** Retain its upper reading guide, percentile graphic and color interpretation. Adjust explanatory wording only as needed to distinguish percentile bars from other units and imported measured results. The guide must not say every number is a percentile or every displayed value comes from DNA once direct measurements are present. This is a copy/layout change, not a scoring change.

Remove **everything on page 3 starting with the heading “Two independent measurements, two reference groups”**, including that explanatory block, the sample-QC strip and the entire **“What this report knew about you, and what it assumed”** block. Replace the released space with a clean **“Explore Your Results”** contents table and the brief navigation copy below.

Relocate rather than discard useful information:

- Known inputs, missing inputs and assumptions move to report section 25.
- The full sample-QC content merges with the existing sequencing-quality/reference material in report section 27, avoiding duplicate tables. Section 25 can link to its compact status summary.
- The measurement/reference explanation moves into the existing methods/reference material in reports 27/29; individual readings retain short context links where useful.

The page-3 contents table links **only to major at-a-glance result sections**. Include the existing summary/highlights and all applicable major overview destinations, omitting a self-link to the reading guide. Do not list individual metrics, long subheadings, technical appendices or every detailed page. Link to each section's first overview page, not directly to its detailed explanation. Generate entries from stable section anchors and final pagination so renamed sections and additional overview pages remain correct.

Use the existing report typography and palette. A balanced two-column table/list is appropriate for the PDF, with ample row spacing, subtle separators, clear section names and right-aligned page numbers. Make the entire row clickable and give each column a clear reading order. On mobile, stack into one column with comfortable tap targets and visible keyboard focus. Keep it on page 3 without shrinking text to fit; use concise labels instead of removing destinations. The remaining part of the report is reached through its readings and the standard navigation controls, rather than an exhaustive contents dump on this page.

Use this short helper text:

**“This report is interactive. Click a section for results at a glance, then a reading for its detailed explanation. Use ‘Back to overview’ to return.”**

Follow it with one short note: **“Detailed explanations and supporting information follow the overview pages.”**

### 12.3 Navigation, integrated displays and completeness

- Every overview reading links to its canonical detail card; every detail card has **Back to overview** targeting its owning section and **Contents** returning to page 3. If a reading is linked from several contexts, the PDF uses its primary home; HTML may also preserve the immediate previous context.
- Page-3 contents is deliberately short. The complete metric index remains available through the report's existing technical/data area and HTML search, with metric ID, label, value, unit, status and destination. Do not turn that index into another top-level results section.
- Per-section counts distinguish unique measurements, repeated contextual views and unavailable analyses. Optional unprovided laboratory fields are not counted as completed microbiome readings.
- Taxonomy search, aliases, organism roles and strain evidence belong within the organism explorer. Condition/context navigation belongs within existing condition sections.
- Last-three-results charts appear beside their relevant reading; full trajectories live with its detail. Evidence/ingredient/formulation comparisons and the printable food checklist live within existing action areas.
- Links from a result's context badge or an applied assumption lead to the exact section-25 row; supplied assay rows link back to the findings they inform. A section-25 record must identify whether it affects a calculation, interpretation, action or context only.
- Methods/units and research evidence remain reachable from relevant readings. Use one canonical value and evidence record wherever it appears.

Every new quantitative graphic shows value, unit, scale meaning, coverage and status. Reference bands use actual reference percentiles or supplied laboratory intervals. Distinguish descriptive high/low, supported favorable/unfavorable interpretation and unavailable/partial evidence using text and icons as well as color. No false healthy middle score for missing input.

The HTML report works on mobile with stacked cards and accessible tables. The PDF has correct bookmarks, named destinations, internal link annotations, repeated table headers and no clipped scientific names or rows. All external HTML links use `target="_blank" rel="noopener noreferrer"`. Rebuild every section/page reference after insertion of section 25; old section 25 becomes 26 and all later section numbers shift consistently. Do not hardcode the old PDF page numbers as final destinations.

### 12.4 One source for displayed values

Render summary, detail, index, HTML, PDF and exports from the same canonical result objects and formatters. Rounding is presentation-only. New aggregate contributor totals reconcile at full precision. Keep unique metrics separate from repeated contexts, genes, organisms, evidence records, recommendations and optional inputs.

Do not make a large count the definition of quality. The complete supported inventory is displayed; unsupported detections are not invented to reach a target. Conversely, a supported result cannot be hidden merely because it has no favorable/unfavorable interpretation.

## 13. CLI, outputs, assets and implementation sequence

### 13.1 Required CLI behavior

Adapt these logical commands to the existing CLI's conventions. Their behaviors and options are required:

```bash
openbiota extend-report --spec 0.8.3 --input results.json --mode cached --out run083
openbiota extend-report --spec 0.8.3 --input results.json --mode incremental --reads-manifest samples.json --out run083
openbiota extend-report --spec 0.8.3 --input results.json --mode full-extension --reads-manifest samples.json --out run083
openbiota extension-assets verify --spec 0.8.3
openbiota extension-reference build --manifest reference_manifest.json --spec 0.8.3
openbiota compare-timepoints --manifest longitudinal_manifest.json --extensions 0.8.3 --out timeline
openbiota import-context --input patient_context.json --sample-id SAMPLE --out normalized_context.json
openbiota import-labs --input laboratory_results.json --sample-id SAMPLE --out normalized_labs.json
openbiota recommend-synbiotics --sample run083/results.ext083.json --mode evidence --out synbiotic_options
openbiota recommend-synbiotics --sample run083/results.ext083.json --mode model-assisted --scenarios scenarios.json --out synbiotic_options
openbiota simulate-substrates --sample run083/results.ext083.json --scenarios scenarios.json --out simulations
openbiota verify-preservation --baseline baseline_results.json --candidate run083/results.json
```

Input manifests are normal application inputs to be implemented, not omitted supplemental build documents. Define their schemas in the repository and generate an example from CLI help. Extend the existing sample manifest/metadata arguments for context intake; the logical import commands may be subcommands or existing import options according to repository conventions. A report with no supplied context or laboratory files must continue to run using declared defaults. Existing commands continue to work unchanged.

Required generated application outputs:

| Output | Contents |
|---|---|
| `results.json` | Existing result objects preserved, with an extension link/object according to backward-compatible schema rules |
| `results.ext083.json` | All new measurements, states, evidence links and companion views |
| `extension_manifest.json` | Input, software, asset, parameter and reference fingerprints |
| `capability_coverage.json` | Every required feature and whether computed, data-limited, externally measured or not applicable |
| `all_metrics.tsv` | All results with IDs, units, statuses and provenance; no hidden truncation |
| `organisms.tsv` | Complete new explorer inventory and taxonomy-resolution state |
| `actions.json` | Deduplicated options, including synbiotic candidates and component comparators, with exact evidence edges and sorting explanations |
| Synbiotic/model objects in existing structured exports | Candidate identities, goal-coverage definitions, ranks/reasons, evidence, contraindication/uncertainty fields and optional model comparators with full provenance |
| Context data in the existing result schema/export | Supplied and effective values, input states, default provenance, laboratory records and exact result/action bindings for report section 25 |
| Report HTML/PDF | Existing content plus new views and navigation |

The user-facing specification remains this one Markdown file. These outputs are artifacts of the implemented application.

### 13.2 Asset manifest and reproducibility

For every downloaded asset store source URL, DOI/publication, accession, version/retrieval date, license/terms, size, SHA-256, file format, expected schema, and purpose. Record corrections/retractions when present. Sources with separate code/data licenses are checked separately. A version string such as `latest` is not a reproducible production manifest.

Core dependencies must remain usable without a restricted commercial knowledge base. Enteropathway, legacy omixer binaries, optional model repositories and some cohort metadata have access/licensing limits; support them through separately eligible adapters rather than making the whole report depend on them. Public raw data and restricted phenotype labels are different assets.

Use bounded downloads, restartable caching and resource estimates derived from the actual chosen databases. Avoid claiming a fixed RAM/disk/runtime budget without measuring the pinned dataset. The cached-only report path must not trigger a multi-terabyte catalogue download.

### 13.3 Implementation order

1. **Freeze preservation fixtures.** Capture all legacy IDs and representative reports.
2. **Build adapters and registry.** Implement extension schemas, manifests, source registry, capability inventory and read-only cache access.
3. **Deliver low-compute additions.** Complete index, taxonomy/target browse, resistance summaries, GABA-breakdown surfacing, new diversity/ratios, metadata/lab import, action aggregation, broad-evidence target/delivery/tradeoff adapters, synbiotic evidence matching/ranking and longitudinal views.
4. **Implement the functional lane.** Build/validate references, run new carbohydrate, acetate/intermediate, vitamin, nitrogen, polyphenol, mucin, bile and gas modules. Add exact P9 analysis and other sequence-resolved host-signaling evidence.
5. **Build new references and research context.** New-method cohort processing, gut-type/aerotolerance/HACK-compatible views, evidence/context tiles and additional capacities. Do not alter legacy references or models.
6. **Integrate the report.** Place new readings in their existing domain overviews/details; rename the functional display headings; rebuild page 3 as the reading guide plus clickable major-section contents; consolidate inputs/assumptions in report section 25; update numbering, bookmarks and all back/detail links.
7. **Validate and release.** Run the preservation suite, analytical fixtures, cohort checks, semantic coverage audit and rendered-report inspection. Optional simulations are clearly isolated and can be enabled independently.

### 13.4 Definition of done

All required additions in A01–A16, including personalized synbiotic evidence matching, the targeted all-evidence recommendation contract in section 19, all 85 indexed TP/BP/TM/BU target-specific seed records (with shared experiments deduplicated), and the detailed coverage appendix are implemented and exercised. Optional model-assisted scenarios in A16 have a working reproducible example and are isolated behind their feature flags. The section-15.7 readiness inventory is resolved for every shipped quantitative capability; no source-identified item is falsely labelled application-validated. Every sequence-computable feature has real code and a tested positive/negative fixture; every data-limited feature has a genuine, specific reason rather than an unimplemented stub. All existing scores pass preservation checks. All new readings have a placement binding to an existing domain section, apart from the combined section 25. Page 3 and every affected section link pass the layout/navigation checks. Every displayed claim, quantitative result and action resolves to its canonical object, implementation method and supporting evidence.

## 14. Acceptance and release tests

Implement all 160 acceptance requirements below as automated checks and the specified rendered-report reviews. Release requires documented passing results, complete feature-to-test traceability, verified data provenance and explicit resolution of failures. A missing implementation, untested reference panel or placeholder output cannot pass as a data limitation. These are acceptance requirements to execute during implementation.

| Test | Required result |
|---|---|
| AT001 | With extension off, canonical legacy results equal the captured baseline, excluding explicitly enumerated volatile run metadata. |
| AT002 | With extension on, existing numerical outputs and labels remain identical for identical baseline inputs. |
| AT003 | Existing disease, GMWI2, donor, strain, pathogen, biofilm, fungal and functional result objects remain unchanged. |
| AT004 | New metadata, reference cohorts or feature values cannot silently enter a legacy model's input graph. |
| AT005 | Existing numerical content is preserved; authorized functional-title changes, page-3 relocation and section-25 insertion match the placement contract. |
| AT006 | Existing CLI entrypoints and consumers of legacy `results.json` remain compatible. |
| AT007 | Cached mode completes without FASTQs/network downloads and lists the precise unavailable extension analyses. |
| AT008 | A new functional database invalidates only dependent extension caches, not unrelated outputs. |
| AT009 | Every required A01–A16 feature and every coverage-appendix metric has an implementation mapping and an exercised fixture; an unimplemented placeholder is not marked data-limited. |
| AT010 | All downloaded assets have actual version, license/source record and verified checksums; no production asset is pinned only to `latest`. |
| AT011 | Fragment multi-mapping and paired-end overlap cannot multiply one molecule's total assignment weight beyond its defined allocation. |
| AT012 | Protein versus nucleotide length normalization and partial-reference handling produce the documented unit, with depth-doubling invariance on simulated data. |
| AT013 | A single generic GH13 hit supports broad amylolysis but does not activate complete resistant-starch utilization. |
| AT014 | Inulin/FOS and lactose/GOS can share evidence without duplicated aggregate fragment counts. |
| AT015 | All 14 carbohydrate views render, with appropriate complete/partial/no-supported-hit/not-assayed states. |
| AT016 | A complete supported alternative B6 route is accepted without requiring the other route; one isolated homolog remains partial. |
| AT017 | B1/B3/B5/B6 synthesis is separated from salvage/uptake and does not generate host vitamin-deficiency text. |
| AT018 | B2/B7/B9/B12/K2 dashboard cells reference unchanged existing result objects. |
| AT019 | AMP-forming `acs` alone does not prove acetate production; assimilation, ADP-forming alternatives and Wood–Ljungdahl CODH/ACS remain distinct. |
| AT020 | GABA production and degradation have separate IDs, values and evidence; a percentile subtraction cannot be labelled net GABA. |
| AT021 | BCAA biosynthesis, amino-acid fermentation and BCFA generation cannot share one biochemical identity. |
| AT022 | Community-distributed route genes do not receive `genome_linked` completeness without linkage evidence. |
| AT023 | UrdA-positive/urolithin-negative fixture does not activate new urolithin or equol modules. |
| AT024 | PQ855390.1-positive fixture supports the specified ucdCFO reaction; a near xanthine-dehydrogenase homolog alone does not. |
| AT025 | Equol component aliases resolve by source/strain identity, not gene-name similarity. |
| AT026 | Glucosinolate core rule accepts BT2158 plus BT2156 or BT2157, with supporting context; generic hydrolases remain nonspecific. |
| AT027 | Mucin enzyme fixtures preserve substrate/linkage specificity and the colonic versus gastric evidence context. |
| AT028 | Hydrogen-producing and hydrogen-consuming classes are separate; high hydrogen potential alone is not an adverse score. |
| AT029 | Bile transformation routes allow supported alternative chemistry; human percentiles cannot be calibrated from rumen cohorts. |
| AT030 | P9 protein/CDS downloads match the specified hashes and coding orientation; Amuc_1831 fails the P9 identity fixture. |
| AT031 | P9 matches are tested against related proteases; unresolved homologs do not become measured secretion or GLP-1 response. |
| AT032 | Historical GBM identifiers retain namespace/version and AND/OR semantics; unresolved mappings cannot silently acquire a current gene identity. |
| AT033 | Tryptamine, serotonin, melatonin and receptor-active derivatives remain distinct chemical entities. |
| AT034 | Generic `adhE` does not yield a high-alcohol strain name, measured ethanol value or auto-brewery diagnosis. |
| AT035 | New reference percentiles are reproducible, participant/study aware and method compatible; raw results remain visible without a reference. |
| AT036 | A sparse zero-valued reference distinguishes whole-cohort prevalence from positive-carrier percentiles. |
| AT037 | For `p=(0.5,0.5)`, new effective Shannon and inverse Simpson both equal 2, and Gini-Simpson equals 0.5. |
| AT038 | For `p=(1)`, effective diversity equals 1; empty profiles return unavailable, not fabricated zeros. |
| AT039 | New genus/family charts are disjoint within each view and reconcile to their declared denominator. |
| AT040 | Zero-denominator ratios are undefined/censored; genus and species ratio definitions are not swapped. |
| AT041 | For `A=.2,N=.6,U=.2`, aerotolerant fraction=.25, trait coverage=.8 and MAPI=`ln(1/3)`; `A=0` yields censored MAPI. |
| AT042 | Missing or conflicting oxygen phenotypes remain unresolved and reduce reported coverage. |
| AT043 | Gut-type assignment uses frozen medoids, feature universe and reference; mixed assignment is allowed and displayed. |
| AT044 | HACK-derived output reproduces the published native definition/test data and distinguishes taxon index from sample score. |
| AT045 | Same-FASTQ reports from different methods cannot be shown as biological improvement over time. |
| AT046 | Identical compositions have Bray–Curtis distance 0; disjoint unit-sum compositions have distance 1. |
| AT047 | Two timepoints yield change, not measured resilience; recovery requires displacement above the declared empirical envelope and is shown only at/after the observed peak. |
| AT048 | Host-filtered reads without provider counts yield unavailable original host fraction; read/pair denominators cannot mix. |
| AT049 | Every taxonomy label in section 16.4 resolves or produces an explicit split/merge/unresolved state; none is silently dropped. |
| AT050 | Species/subspecies and alias identities are not double counted; no fixed species-count target creates false detections. |
| AT051 | Named-target coverage includes all targets in section 7.2 with actual installed/reference/assay/execution states. |
| AT052 | DNA-only input cannot report an assessed negative for the listed RNA viruses. |
| AT053 | A species detection does not automatically receive an unobserved toxin, pathotype, commercial probiotic strain or resistance phenotype. |
| AT054 | AMR class totals deduplicate alleles, family aliases and parent/child drug categories; ambiguous carriers remain explicit. |
| AT055 | New action cards retain animal/ex-vivo/in-vitro evidence visibly; lack of an RCT alone never removes a supported research option. |
| AT056 | Combination, formulation, strain and host-species differences prevent unsupported transfer of an intervention claim. |
| AT057 | Null/contrary study findings are retained and cannot be converted into positive target-suppression evidence. |
| AT058 | A study concentration or animal dose cannot become a human regimen through automatic unit conversion. |
| AT059 | Reported allergy/intolerance changes new action prioritization; unknown recorded facts remain distinct from visible effective defaults and are never relabelled as patient-confirmed negatives. |
| AT060 | Duplicate food/action items merge with all supported targets and tradeoffs; overlapping targets do not create a clinical net-benefit score. |
| AT061 | Imported laboratory values retain identity, method, matrix, date, censoring and original units; supported inputs link to relevant findings/actions and section-25 provenance; DNA proxies cannot occupy measured-analyte fields. |
| AT062 | With no external labs, section 25 displays “Not supplied”; the report runs with declared context defaults and contains no fabricated normal/zero assay values or repeated empty laboratory charts. |
| AT063 | Model scenarios report coverage, assumptions, feasible status and uncertainty, and never modify protected scores. |
| AT064 | Missing food-DNA evidence cannot be interpreted as zero dietary intake. |
| AT065 | Each new displayed value matches its canonical object across summary, detail, index, HTML, PDF and export. |
| AT066 | Existing domain overviews and the complete deeper index cover every supported reading; page-3 contents lists only major overview destinations and is not an exhaustive metric index. |
| AT067 | Mobile cards, long taxonomy names, wide evidence tables and printable action lists render without clipping. |
| AT068 | PDF bookmarks and internal links reach correct pages; external HTML links have the required target/rel attributes. |
| AT069 | New numerical graphics distinguish scale meaning, reference position, evidence coverage and health interpretation without colour alone. |
| AT070 | All new claims/actions resolve to embedded source IDs/URLs; source corrections and access restrictions are preserved. |
| AT071 | Every A01–A16 output has an explicit existing-section placement or the authorized section-25 context placement; no standalone feature chapters are created from capability IDs. |
| AT072 | Cover and report sections 10/20 consistently use the new Microbial Functions & Metabolism umbrella; all previous and new functions remain discoverable with unchanged canonical values. |
| AT073 | Page 3 contains none of the old content beginning with “Two independent measurements, two reference groups”; that reference explanation, QC and assumptions are present at their specified new destinations. |
| AT074 | The complete major-overview contents fits page 3 with readable typography; every row reaches its actual overview in both HTML and PDF, and no detailed-page inventory is included in that table. |
| AT075 | Every new reading/detail pair has working detail, Back to overview and Contents links; PDF destinations and HTML keyboard/touch behavior are verified. |
| AT076 | Report section 25 is Your Information & Report Context; old section 25 becomes 26, old26→27, old27→28 and old28→29, with all references/bookmarks updated. |
| AT077 | All influential supplied inputs and assumptions appear in section 25 with effective value, provenance, rationale and consuming result/action; a missing optional input does not block unrelated analyses. |
| AT078 | A supplied value replaces its applicable default through declared dependencies; the recorded fact and prior assumption remain distinguishable, and prior reports remain immutable. |
| AT079 | The baseline existing-default behavior is preserved; new defaults are versioned and testable, and neither invented diagnoses nor invented numerical laboratory values are used as fallbacks. |
| AT080 | A compatible supplied SCFA or other mapped assay appears beside the associated functional result as a separately labelled measurement and links to its section-25 record; no new standalone laboratory chapter appears. |
| AT081 | Wrong-participant, incompatible-matrix, conflicting or stale imports cannot silently change another sample's interpretation; their status and any allowed use are explicit. |
| AT082 | Context-only versus calculation/interpretation/action use is visible; a record without an implemented consumer is labelled recorded/not used rather than falsely integrated. |
| AT083 | Time trends, scenario views, action summaries, strain details and AMR additions stay within their existing domain sections rather than generating a second overview or new top-level feature sections. |
| AT084 | New layout regressions are checked on the report with no optional context and with a populated context/lab fixture, including page 3, the expanded functional overviews, section 25, mobile layout and final pagination. |
| AT085 | Core synbiotic evidence mode works from valid cached findings without a metabolic solver, private cohort or commercial product API. |
| AT086 | All SYN001–SYN017 records are ingested or carry a precise incomplete-field state; sources, identities, preparation, comparator and outcome type are preserved. |
| AT087 | IVS-1 and BB-12 GOS factorial fixtures retain no demonstrated added combination benefit; component benefits cannot flip the interaction result positive. |
| AT088 | Bi-07/XOS preserves its component-only comparators and absence of added strain enrichment; missing arms in other trials remain missing. |
| AT089 | PLMB0001 adult evidence cannot map to PBI001 mouse evidence; LAHUC cannot map to DSM14662/L1-92. |
| AT090 | A species-only sample detection, a supplier label and high genome ANI cannot independently resolve a commercial probiotic strain. |
| AT091 | The DS-01 trial is one 24-strain bundle record with AFU units; it cannot generate 24 independently proven pair effects or a personalized allocation claim. |
| AT092 | Measured growth-table observations and predicted 3,083-genome utilization entries have distinct evidence types; unknown and tested-negative remain different. |
| AT093 | Substrate aliases preserve chain length, mixtures, residual sugars and hydrolysis; intact beta-glucan cannot silently inherit hydrolysate findings. |
| AT094 | No-addition, substrate-only and probiotic-only alternatives are retained; fibre-only is allowed to outrank every combination. |
| AT095 | A reported allergy or explicit contraindication changes start-here eligibility with a cited reason; the research record remains visible. |
| AT096 | Animal, infant, ex-vivo and adult outcomes retain population/model applicability and cannot silently become the same human clinical evidence. |
| AT097 | With equal weights and two of three supported goals, goal coverage is 66.666…%; one pair-supported goal is 33.333…%; neither field is labelled efficacy probability. |
| AT098 | No goals or all-zero goal weights yield null coverage; negative/nonfinite weights are rejected; contrary evidence is visible and unknown tradeoffs cannot masquerade as measured zero harm. |
| AT099 | Ranking is deterministic for identical inputs/policy; preference and goal changes recompute only dependent new actions and record their reasons. |
| AT100 | Novel multi-component candidates do not inherit whole-bundle efficacy or sum pairwise effects; dropping components and single-option comparisons are supported. |
| AT101 | Probiotic input perturbations preserve total abundance one; .95 plus five .04 additions fails the mass fixture rather than entering production silently. |
| AT102 | Model mapping reports retained abundance and threshold losses on the declared bacterial denominator; unmapped organisms are not reassigned to convenient models. |
| AT103 | Missing target exchange or nonoptimal solver status returns unavailable with its reason, not zero flux or a successful prediction. |
| AT104 | Four-arm values 1,3,4,8 produce pair delta7, over-substrate4, over-probiotic5 and interaction2; missing any arm prevents that contrast. |
| AT105 | Scenario media record added-versus-substituted carbon, exchange IDs and units; model uptake caps cannot become oral doses. |
| AT106 | Strict OOXML S5 ingestion succeeds or emits a specific format failure; it cannot accept zero sheets as an empty valid dataset. |
| AT107 | The section-11.7 S5 sample reproduces the raw-flux fixture and its fibre-only winner with declared numerical tolerance. |
| AT108 | The same sample has maltodextrin first under the separate growth-normalized objective; raw and normalized flux units/order cannot be swapped. |
| AT109 | A published-output replay and an actually solved model example have different execution states and provenance; neither is called a measured user outcome. |
| AT110 | The 154-input/156-output discrepancy is audited by sample ID; missing inputs and conflicting growth thresholds are not silently repaired with fabricated data. |
| AT111 | Predeclared scenario sensitivity displays a range/rank stability, not an empirical confidence interval or treatment-success probability. |
| AT112 | The MiRNN adapter accepts its defined culture domain and cannot present arbitrary-patient metagenome input as validated personalization. |
| AT113 | Strain/participant/study-aware validation keeps synonyms, repeated participants and technical wells out of both training and held-out sets. |
| AT114 | Synbiotic options render inside report15 and mechanism details inside report20 with working links; no new top-level section is created. |
| AT115 | Every A01–A16 feature exposes separate source/retrieval/registry/reproduction/application readiness; not-implemented is not labelled scientifically unresolvable. |
| AT116 | New percentile/model outputs require their actual frozen reference/trained artifact; missing references do not suppress supported raw measurements. |
| AT117 | F17 protein retrieval matches 624 aa and its normalized sequence hash; the candidate enzyme is not labelled serotonin-specific based only on tryptamine activity. |
| AT118 | F17 reference MAG and reconstructed study construct remain distinct; a retrieved homolog cannot become a verified WL68 genome. |
| AT119 | The missing ABS author FASTA/metadata are explicitly tracked; an independent ethanol panel cannot claim exact author-panel reproduction. |
| AT120 | Source file checksums, conditional access, archive-specific terms and mutable URLs are recorded; an HTTP-success check alone cannot mark an asset installed/validated. |
| AT121 | A target-matched controlled culture experiment produces a normal finding-specific candidate without any human study; default report/API retrieval contains it. |
| AT122 | Animal, biochemical, ex-vivo, organoid, culture and human fixtures preserve their setting and actual endpoint; no setting is relabelled clinical efficacy. |
| AT123 | An unrelated human trial cannot promote or replace a directly target-matched inhibition observation; ranking reasons expose the target/endpoint joins. |
| AT124 | A laboratory-supported compound with unknown intestinal exposure remains visible as an option with unresolved exposure, rather than “ineffective” or contraindicated. |
| AT125 | An enteric-delivery hypothesis is retained, but cannot become a verified retail preparation or measured colonic concentration without its own evidence. |
| AT126 | Pure allicin, standardized garlic preparations, aged garlic extract and garlic oil remain distinct; neither efficacy nor bioavailability is transferred by ingredient-name similarity alone. |
| AT127 | Growth inhibition, bacterial killing, prevention of attachment, mature-biofilm biomass, viable-cell reduction and dispersal have separate endpoint fixtures; crystal-violet reduction cannot become eradication. |
| AT128 | Existing biofilm scores retrieve related options without claiming physical biofilm quantity was measured or changing the score. |
| AT129 | Existing high/low/undetected organism, pathogen, mycobiome, function and symptom findings each have a typed option-lookup record and a working action link. |
| AT130 | “Not searched,” failed search, incomplete extraction, no matching record, conflicting evidence and supported options are distinct states; none automatically renders “untreatable.” |
| AT131 | The complete candidate count and expandable list agree with canonical output; the top-three view does not discard additional or preclinical options. |
| AT132 | AOR/TO-A and CBM588/MIYAIRI588 records cannot exchange strain-specific efficacy or survival evidence. |
| AT133 | A partial 16S accession cannot pass a whole-genome identity check; a manufacturer sequence link is verified for accession type. |
| AT134 | Verified T-110/TO-A/TO-A product identities, constituent quantities and lactose carrier remain attached to the correct market-label snapshot; single-component evidence is not a whole-bundle result. |
| AT135 | Live recovery during ingestion, sequence detection, washout persistence and human benefit are independent fields; absent permanent engraftment cannot eliminate a supported option. |
| AT136 | The exact TO-A culture butyrate edge and animal consortium butyrate edge retain different preparations, hosts and endpoints. |
| AT137 | Strain-specific HMO growth positives and tested negatives are preserved; 2′FL/3FL/3′SL/LNT/LNnT cannot collapse to one “HMO selectivity” flag. |
| AT138 | A branded HMO product with undisclosed constituent structures cannot inherit a specific purified-HMO experiment; composition uncertainty is displayed without discarding relevant category-level research. |
| AT139 | R. gnavus monoculture inability to consume intact starch cannot erase measured cross-feeding via another organism's starch breakdown products. |
| AT140 | Glutamine-related genes or generic nitrogen requirements cannot independently label oral glutamine as feeding an undesirable organism or contraindicate it. |
| AT141 | Beneficial and unfavorable ecological effects of the same ingredient remain separate context-specific edges, including conflicting berberine directions where documented. |
| AT142 | KPV cell/animal observations and BPC-157 preclinical observations remain visible as investigational candidates with actual delivery route; no human efficacy, safety or dosing is fabricated. |
| AT143 | Zinc-carnosine and glutamine host-permeability observations do not relabel stool sequence proxies as a direct permeability measurement. |
| AT144 | Exact-strain competitive-colonization prevention cannot become clearance of established carriage; challenge timing and antibiotic conditioning are retained. |
| AT145 | Positive and null trials for the same exact probiotic/target remain independently retrievable; a favorable-only import fails validation. |
| AT146 | NPASS/ChEMBL/PubChem copies of the same experimental observation count once; independent replication excludes shared-cohort publications and technical replicate wells. |
| AT147 | ChEMBL whole-organism target-assignment confidence 1 is not rejected by a protein-target confidence cutoff; database confidence is never displayed as therapeutic efficacy. |
| AT148 | MIC values with inequality signs, units, salt forms and media remain intact; no conversion to a human oral dose occurs. |
| AT149 | Docking, genome prediction and simulated metabolism are visible as mechanisms/models, never experimentally observed antimicrobial activity. |
| AT150 | A supplied allergy/interaction can exclude a candidate from start-here with an explicit reason; unknown patient facts and absent trials cannot do so. |
| AT151 | Compound-first and target-first discovery paths operate across all evidence settings; no implicit human-only/RCT filter is present in default queries. |
| AT152 | Cached recommendation rendering succeeds offline; source timeouts, changed schemas and access failures produce retrieval errors, not negative evidence. |
| AT153 | Every source-derived card resolves to an exact publication/assay locator and supported extraction level; missing concentrations or strain identifiers stay null. |
| AT154 | Unknown product availability and experimental delivery are displayed separately from effect evidence; ingredient-level options work without affiliate/retail links. |
| AT155 | The version-2 ranker is deterministic; unrelated targets cannot gain rank via favorable studies, and incomplete tradeoff assessment cannot be treated as zero harm. |
| AT156 | Unsupported simultaneous combinations do not inherit additive efficacy or synergy; component antagonism and unresolved interactions remain visible. |
| AT157 | All embedded new resource/seed rows import with unique local IDs, provenance and explicit unresolved fields; no research-note path is required at runtime or for implementation. |
| AT158 | An obsolete/unmaintained source is dated and labelled; a guide whose official site is on hold cannot be presented as a current authoritative complete catalogue. |
| AT159 | Existing metric values and protected computation fingerprints remain unchanged; recommendation/ranking changes are confined to the authorized new planner and its evidence records. |
| AT160 | Rendered finding cards show target, preparation, evidence setting, actual effect, delivery and tradeoffs clearly, with all primary links and existing section navigation functioning. |

Analytical validation must include held-out positive and near-negative reference genomes, synthetic reads matched to real read length and depth, mixtures, dilution series and relevant public mocks. Report precision/recall and the actual evaluable denominator. Tune new-marker thresholds on development data and test them on held-out strains/cohorts; do not calibrate them to make the supplied samples look better. Research estimates remain useful when their uncertainty is visible.

## 15. Embedded research, data and code register

The source groups below are implementation inputs. Retrieve the cited primary materials, preserve the exact release and terms, and bind each imported sequence/definition to its source. Asset access limitations are stated where identified. This register and the intervention matrix are part of this specification.

### 15.1 Functional sources

**F01** — dbCAN3 original methods: https://pmc.ncbi.nlm.nih.gov/articles/PMC10320055/; DOI 10.1093/nar/gkad328.

**F02** — dbCAN-seq CGC/PUL substrate expansion: https://pmc.ncbi.nlm.nih.gov/articles/PMC9825555/.

**F03** — official V5 code: https://github.com/bcb-unl/run_dbcan; GPL-3.0 code.

**F04** — exact database documentation: https://run-dbcan.readthedocs.io/en/latest/getting_started/database_description.html; substrate rules https://run-dbcan.readthedocs.io/en/latest/user_guide/predict_CGC_substrate.html .

**F05** — open data distribution: https://registry.opendata.aws/run_dbcan/; declares no restrictions on use of data; https://pro.unl.edu/dbCAN2/ latest HMM release. Optional SignalP6/DeepTMHMM binaries are not bundled and have separate installation/license terms.

**F06** — Magnúsdóttir et al., Systematic genome assessment of B-vitamin biosynthesis (2015): https://pmc.ncbi.nlm.nih.gov/articles/PMC4403557/; DOI 10.3389/fgene.2015.00148. Curated eight-vitamin reconstructions across 256 gut bacteria.

**F07** — Rodionov et al., Micronutrient Requirements and Sharing Capabilities of the Human Gut Microbiome (2019): https://pmc.ncbi.nlm.nih.gov/articles/PMC6593275/. Curated biosynthesis/salvage/uptake for B1/B2/B3/B5/B6/B7/B9/B12 and queuosine over 2,228 genomes, 690 cultured species; supplementary pathway/gene matrices are the import source. These genome phenotype priors support validation; read evidence still determines individual sample output.

**F08** — Vital et al. butyrate pathways: https://pmc.ncbi.nlm.nih.gov/articles/PMC3994512/; DOI 10.1128/mBio.00889-14. Four alternative upstream pathways (acetyl-CoA, glutarate, lysine, 4-aminobutyrate) and different terminal enzymes. Use published sequences and functional assignments, not `but` alone for all capacity.

**F09** — Reichardt et al. propionate pathways: https://pmc.ncbi.nlm.nih.gov/articles/PMC4030238/; DOI 10.1038/ismej.2014.14. Succinate, acrylate and propanediol routes; diagnostic candidates mmdA, lcdA, pduP need pathway context and near-homolog discrimination. Existing propionate output remains untouched; new pathway breakdown makes the chemistry clear.

**F10** — Dodd et al. aromatic routes (2017): https://www.nature.com/articles/nature24661. Genetics establishes a shared aromatic reductive pathway; FldABC/FldH/FldI/AcdA evidence must not be called IPA-exclusive without substrate context.

**F11** — Liu et al. Stickland metabolism (2022): https://www.nature.com/articles/s41564-022-01109-9; PMC9089323. Pathway and knockout data distinguish fermentation from de novo biosynthesis.

**F12** — Phenylacetate/PAGln microbial versus host steps: https://pmc.ncbi.nlm.nih.gov/articles/PMC9839529/. PorA/PPFOR pathways produce microbial precursors; host conjugation creates PAGln. Report precursor capacity, not circulating PAGln.

**F13** — p-cresol hpdBCA experimentally tested: https://pmc.ncbi.nlm.nih.gov/articles/PMC7925072/; DOI 10.1128/JB.00282-20. Substrate-driven expression important. New 2026 upstream route: https://pmc.ncbi.nlm.nih.gov/articles/PMC13285125/ adds tyrosine→HPA biochemistry; preserve HpdBCA versus generic decarboxylases.

**F14** — Multiple substrate sources of TMA: https://pmc.ncbi.nlm.nih.gov/articles/PMC8364193/; https://www.nature.com/articles/s41564-021-01010-x . Add anaerobic gamma-butyrobetaine `bbu`/`gbu` route details as aliases with validated exact accession mappings, and split oxygen-dependent CntAB from anaerobic route. TMAO is host-converted; no stool gene measure of serum TMAO.

**F15** — Valles-Colomer et al. 2019 GBMs: https://www.nature.com/articles/s41564-018-0337-x; https://raeslab.org/software/gbms.html .

**F16** — https://github.com/omixer/omixer-rpmR , commit `184f1cc99aefed722e20eb00eda082348a064c4e`; files `inst/extdata/GBMs.v1.0.names` and `GBMs.v1.0.txt`. README wrapper GPL3, **bundled omixer-rpm.jar academic noncommercial** . Do not redistribute under blanket open commercial claim.

**F17** — Moretti 2025, DOI https://doi.org/10.1016/j.celrep.2025.116434 ; full primary text https://pub.epsilon.slu.se/38878/1/moretti-c-h-et-al-20251117.pdf . The original Ls consortium supports 5-HTP→serotonin; isolated *L. ruminis* WL43/*L. mucosae* WL68 and reconstituted cultures are not interchangeable with it. Raw project **PRJEB96399** has 211 public FASTQ run records at source verification. Candidate protein **WP_407418999.1** is currently retrievable despite the paper's removed-record warning: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=WP_407418999.1&rettype=fasta&retmode=text . Its 624-aa uppercase unwrapped sequence SHA256 is `3b161781132370598a9462d8d0573589ca4c72d5ee04bbb8c88cb3bab9ffcf62`. Reference contig **NZ_JAZOPK010000041.1**, assembly **GCF_045655505.1**, BioSample **SAMN37894531**, belongs to a reference MAG, not a deposited WL68 isolate. The reconstructed study CDS was 1,875 bp with a 60-bp reference-derived gap fill; the exact study construct remains to be resolved. Purified-enzyme evidence uses tryptophan and produces tryptamine; serotonin-specific activity is not established by that assay. Preserve candidate-enzyme versus consortium evidence. Human clinical metadata are request/approval-dependent. Do not invent a validated serotonin marker.

**F18** — Tingler et al. (online Dec 2025 / issue Feb 2026): https://pmc.ncbi.nlm.nih.gov/articles/PMC12860712/; DOI 10.1016/j.isci.2025.114424. In vitro eight-commensal panel produced diverse compounds, but none made serotonin or norepinephrine under tested conditions. Context and strain specificity matter; broad microbial neurotransmitter marketing is not a calibrated test.

**F19** — Pidgeon et al. (2025): https://www.nature.com/articles/s41467-025-56266-2; experimentally validated ucdCFO and public operon sequence. DNA prevalent but expression determines ex vivo conversion.

**F20** — additional independent 2025 urolithin enzymes: https://www.pnas.org/doi/10.1073/pnas.2501312122; novel regioselective pathways https://pmc.ncbi.nlm.nih.gov/articles/PMC12658268/. Import only validated sequences and distinguish molecular products.

**F21** — equol purified enzymes: https://pubmed.ncbi.nlm.nih.gov/22286043/; https://journals.asm.org/doi/10.1128/aem.00410-12; Transcript/cluster support https://pmc.ncbi.nlm.nih.gov/articles/PMC6566806/.

**F22** — UrdA near-homolog discrimination: https://www.nature.com/articles/s41467-021-21548-y; https://pmc.ncbi.nlm.nih.gov/articles/PMC7676231/. Align reference numbering: published position 373 Y/M vs H distinguishes supported UrdA from fumarate reductase homologs; a raw read offset 373 is invalid. Incomplete coverage of distinguishing residues means unresolved specificity.

**F23** — Liou et al. (2020), A Metabolic Pathway for Activation of Dietary Glucosinolates by a Human Gut Symbiont: https://pubmed.ncbi.nlm.nih.gov/32084341/; primary https://www.sciencedirect.com/science/article/pii/S0092867420301045. Genetic transfer, knockout, enzyme and mouse assays.

**F24** — Welsh et al. (2025): https://www.nature.com/articles/s41564-025-02154-w; Group B [FeFe]-hydrogenase dominates healthy gut H2 production; study uses HydDB, DIAMOND, identity/length filters and 14 single-copy ribosomal markers for genome-equivalent normalization. Data **PRJEB45397**, **PRJEB70412**. **High H2 potential is not automatically bad.**

**F25** — HydDB current official repository https://github.com/GreeningLab/HydDB; old Aarhus site unmaintained. Pin HMMs, sequence classifier and citation separately; missing explicit license in repository landing page requires verify before bundling, not before using published biochemical facts.

**F26** — taurine/isethionate sulfide route: https://pmc.ncbi.nlm.nih.gov/articles/PMC6386719/ , IslA plus cognate activase IslB and sulfite reduction, not generic GRE or Bilophila abundance. Sulfoacetate route SauCD/SauS/TauF: https://pmc.ncbi.nlm.nih.gov/articles/PMC10413351/. B. wadsworthia is sulfite-reducing; do not call it a sulfate reducer from genus identity.

**F27** — cysteine sulfide mechanisms: https://pubmed.ncbi.nlm.nih.gov/34745023/; curated primary/secondary/erroneous enzymes across UHGG. Gene specificity and culture context must be retained. Stool potential cannot replace breath tests or diagnose SIBO/intestinal methanogen overgrowth.

**F28** — Liu et al. (2023): https://pmc.ncbi.nlm.nih.gov/articles/PMC10421625/; DOI 10.1016/j.cell.2023.06.010; data **PRJNA850627**, paper Data S1 gut contigs.

**F29** — Kasahara et al. (2023): https://pmc.ncbi.nlm.nih.gov/articles/PMC10311284/; data **PRJNA904303**, **PRJNA911264**, **PRJNA903666**, controlled human **EGAS00001004480**.

**F30** — CAPA spermidine route: https://pmc.ncbi.nlm.nih.gov/articles/PMC10599626/; DOI 10.1126/sciadv.adj9075.

**F31** — 2025 polyamine stable-isotope/multi-omics: https://pmc.ncbi.nlm.nih.gov/articles/PMC11812404/; public validation can use IBDMDB metatranscriptome/metabolome.

**F32** — 2026 E. rectale glutathione: https://link.springer.com/article/10.1186/s40168-026-02457-y; canonical fusion enzyme biochemical primary https://pmc.ncbi.nlm.nih.gov/articles/PMC1112035/.

**F33** — bai pathway biochemical reconstitution: https://www.nature.com/articles/s41586-020-2396-4; (six required catalytic enzymes for CA→DCA; transporter/context tracked separately).

**F34** — BSH is dual hydrolase/acyltransferase, two independent 2024 Nature papers: https://www.nature.com/articles/s41586-023-06990-w and https://www.nature.com/articles/s41586-024-07017-8;

**F35** — 2026 direct conjugated-BA transformation: https://www.nature.com/articles/s41467-026-68556-4; published 2026-02-23. Some HSDHs transform conjugated BAs before deconjugation; BSH is not prerequisite for every possible secondary transformation.

**F36** — **BileActome** (2026, 27 HMM families) https://github.com/1362996609zby/BileActome; DOI 10.1093/ismeco/ycag205; commit `2641f51eed97c6880adaedf87ab374033f84850c` (2026-08-16), CC BY4.0 original repo material. Download `database file/*.hmm` and example scripts. Families: BSH_g/BSH_t; 3α/β,7α/β,12α/β HSDH subtypes; BaiA1/A2/B/CD/E/F/G/H/I/J/K/N/O/P. HMM default full-sequence E≤1e-3 is the authors' starting rule, not a clinical validation threshold. **Rumen validation, not a human healthy-reference cohort**: validate human sequences using F33–F35; do not import rumen abundance norms.

**F37** — HUMAnN official https://github.com/biobakery/humann. Supports stratified gene families, reactions, pathways and pathway coverage with community/taxon contributions. Pin compatible tool and database versions; do not simply mix HUMAnN3 ChocoPhlAn with arbitrary MP4 output without documented compatibility. Existing tables are reused only when provenance matches.

**F38** — Enteropathway (2024): https://pmc.ncbi.nlm.nih.gov/articles/PMC11367760/; DOI 10.1093/bib/bbae419; 3,269 compounds, 3,677 reactions, 876 modules. **Not unconditional commercial open data:** current 2026-07-07 official terms https://comp.life.isct.ac.jp/enteropathway_info/ CC BY-NC-SA4.0, commercial use paid agreement, academic-login download/API. Use as optional properly licensed adapter/research navigator, not mandatory bundled database. Core independently curated routes from public papers must remain usable without it.

**F39** — GUS substrate classes: https://pmc.ncbi.nlm.nih.gov/articles/PMC5533298/ and https://pmc.ncbi.nlm.nih.gov/articles/PMC6351562/; Add substrate/loop-structure detail to a new beta-glucuronidase table; total GUS genes cannot specify estrogen reactivation. More recent atlas https://www.nature.com/articles/s41467-025-65679-y requires GUS domains plus catalytic/recognition motifs, not generic GH2 alone.

**F40** — Digoxin mechanism source clarification for new independent drug-metabolism evidence card: https://pmc.ncbi.nlm.nih.gov/articles/PMC3736355/ (Haiser2013) and https://pmc.ncbi.nlm.nih.gov/articles/PMC5953540/ (Koppel2018). **Arginine inhibits**, not increases, Cgr-mediated digoxin reduction in these experiments. Use as biochemical context, not an individual drug concentration or dose instruction.

**F41** — Bakshani et al. 2025: https://www.nature.com/articles/s41564-024-01911-7; https://pmc.ncbi.nlm.nih.gov/articles/PMC11790493/;

**F42** — Dey et al., 27 July 2026: https://www.nature.com/articles/s41564-026-02424-1. Colonic-mucin-specific sulfatase work, including Amuc1755 and Amuc0953, cellular localization, mutations and structural validation. Eight of ten S1 sulfatases were kinetically characterized; substrate presentation matters. Use this study's corrected version (publisher correction 9 September 2026, https://doi.org/10.1038/s41564-026-02504-2). A. muciniphila mucin utilization alone supports neither blanket favorable nor unfavorable direction. Sequence potential cannot establish excessive host mucus consumption.

**F43** — GLP-1 mechanistic primary sources: Tolhurst et al. 2012, https://pubmed.ncbi.nlm.nih.gov/22190648/ and https://pmc.ncbi.nlm.nih.gov/articles/PMC3266401/ (, ), SCFA-triggered host GLP-1 secretion via FFAR2 in experimental systems; Chimerel et al. 2014, https://pubmed.ncbi.nlm.nih.gov/25456122/ and https://pmc.ncbi.nlm.nih.gov/articles/PMC4308618/ (, ), indole effects differ with exposure duration, including prolonged reduction in GLP-1 secretion. Thus even an 'indole activator' one-direction label would be inaccurate. Link each source in the new host-signaling explanation rather than manufacturing a hormone concentration from microbial genes.

### 15.2 Ecology, reference and interpretation sources

| ID | Research/resource and implementation role | Primary links and available assets |
|---|---|---|
| E01 | Arumugam 2011, enterotypes; foundational community-configuration analysis | https://www.nature.com/articles/nature09944 |
| E02 | Costea 2018, enterotypes across the composition landscape; supports mixed/continuous interpretation | https://www.nature.com/articles/s41564-017-0072-8 |
| E03 | Vandeputte 2017, quantitative microbiome profiling; relative abundance versus microbial load | https://www.nature.com/articles/nature24460 |
| E04 | 2018 gut-redox/MAPI work; aerotolerant-to-anaerobe log ratio | https://www.sciencedirect.com/science/article/pii/S2452231718300150 |
| E05 | 2025 MAPI assessment; research ecological association, not direct oxygen measurement | https://doi.org/10.1128/spectrum.00271-25 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC12323360/ |
| E06 | BacDive 2025; cultured strain phenotypes, including oxygen tolerance | https://academic.oup.com/nar/article/53/D1/D748/7848838 ; https://api.bacdive.dsmz.de/ ; https://hub.dsmz.de/wiki/bacdive/api/ |
| E07 | metaTraits 2026 and porTraits; large genome/trait resource, with measured versus predicted evidence retained | https://academic.oup.com/nar/article/54/D1/D835/8343513 ; https://metatraits.embl.de/ ; https://github.com/grp-bork/porTraits |
| E08 | curatedMetagenomicData 3; cohort discovery, participant metadata, compatible profiles and new-reference construction | https://www.nature.com/articles/s41467-025-66888-1 ; https://github.com/waldronlab/curatedMetagenomicData ; https://waldronlab.io/curatedMetagenomicData/articles/curatedMetagenomicData.html |
| E09 | Asnicar 2025 online/2026 issue, gut microorganisms associated with health/nutrition; published 661-SGB association ranks | https://www.nature.com/articles/s41586-025-09854-7 ; Supplementary Table 5; public profiles/basic metadata https://doi.org/10.5281/zenodo.15307999 ; analysis archive https://doi.org/10.5281/zenodo.17236261 ; https://github.com/SegataLab/inverse_var_weight ; raw projects PRJEB75460, PRJEB75462, PRJEB75463, PRJEB75464 |
| E10 | HACK 2025; health/core/stability-associated taxon index and distinct top-17 sample score | https://doi.org/10.1016/j.celrep.2025.115378 ; https://github.com/tsg-microbiome/HACK_index_2024 ; https://doi.org/10.5281/zenodo.14742683 |
| E11 | Mazzoni 2026, host-derived fecal DNA/inflammation; use corrected article and prefilter counts | https://doi.org/10.1186/s40168-026-02344-6 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC13045140/ ; microbial reads PRJNA1265906 are human-filtered and cannot reproduce the original host fraction; some adult data are controlled |
| E12 | BactoTraits 2026; additional curated bacterial trait evidence | https://www.nature.com/articles/s41597-026-06652-2 |
| E13 | Costea 2017, processing standards; specimen/method comparability | https://doi.org/10.1038/nbt.3960 |
| E14 | Daily quantitative microbiome variability and stool-marker variability; repeated-sample context | https://www.nature.com/articles/s41467-021-27098-7 ; https://www.nature.com/articles/s41598-024-75477-z |
| E15 | Lipid-A/LPS structure and context-specific immune effects; validate new structure/interpretation fields | https://pmc.ncbi.nlm.nih.gov/articles/PMC4950857/ ; https://doi.org/10.1128/mSystems.00046-17 ; https://doi.org/10.1038/s41564-025-01930-y |
| E16 | MEDI final 2025 food-DNA study; optional dietary-exposure context | https://www.nature.com/articles/s42255-025-01220-1 ; https://github.com/Gibbons-Lab/medi ; https://github.com/Gibbons-Lab/medi-paper/blob/main/db/data/manifest.csv |
| E17 | Prevotella:Bacteroides diet-response studies; descriptive ratio context, not universal cutoffs | https://www.nature.com/articles/ijo2017220 ; https://www.nature.com/articles/s41366-018-0093-2 |
| E18 | Species-specific colorectal qPCR ratio research; do not substitute genus-level shotgun ratios | https://doi.org/10.1373/clinchem.2018.289728 ; https://pubmed.ncbi.nlm.nih.gov/29914865/ |

**E09 import details:** native released profiles use a Jan21 SGB database and beta MetaPhlAn code. Preserve the source taxonomy and resolve its relationship to the current catalogue. Public profiles/basic metadata do not include unrestricted access to all clinical endpoints. Use the published rank manifest for an association annotation, not a newly claimed clinical reference population.

**E10 import details:** the species-level HACK index and the top-17 sample score are different artifacts. Fetch the exact species manifest and native ranking algorithm from the archive, reproduce the released examples, then define a frozen external reference for new samples. Do not rank a sample against only the five in-house samples or a changing batch of donors. Any adaptation is labelled HACK-derived with its own version and validation.

### 15.3 Context and host-signaling sources

| ID | Primary evidence and implementation role | Links/assets |
|---|---|---|
| C01 | Chambers 2015, targeted colonic propionate and human appetite/PYY/GLP-1 responses | https://doi.org/10.1136/gutjnl-2014-307913 ; https://pubmed.ncbi.nlm.nih.gov/25500202/ |
| C02 | Yoon 2021 P9 experiments; Di 2024 heterologous P9 expression; independent genomic mapping | https://www.nature.com/articles/s41564-021-00880-5 ; https://doi.org/10.1007/s11274-024-04012-z ; https://doi.org/10.1186/s12866-021-02415-8 ; original sequence-identifying supplement https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41564-021-00880-5/MediaObjects/41564_2021_880_MOESM3_ESM.xlsx |
| C03 | Hsu 2026, gut microbial ethanol metabolism in an auto-brewery observational cohort | https://www.nature.com/articles/s41564-025-02225-y ; https://pmc.ncbi.nlm.nih.gov/articles/PMC13244660/ ; PRJNA1257465 ; https://github.com/hlfreund/ABS_CHsu |
| C04 | Yuan 2019, high-alcohol-producing Klebsiella strains and fatty-liver experiments | https://doi.org/10.1016/j.cmet.2019.08.018 ; https://pubmed.ncbi.nlm.nih.gov/31543403/ |
| C05 | Meijnikman 2022, microbial ethanol and human fatty-liver context | https://www.nature.com/articles/s41591-022-02016-6 |
| C06 | 2025 microbial metabolites, amygdala excitability and anxiety-linked behavior; experimental mechanism | https://doi.org/10.1038/s44321-024-00179-y ; https://pmc.ncbi.nlm.nih.gov/articles/PMC11821874/ |
| C07 | Su 2024, pediatric multikingdom/functional autism markers; corrected source and metadata access | https://www.nature.com/articles/s41564-024-01739-1 ; corrections https://www.nature.com/articles/s41564-024-01900-w and https://www.nature.com/articles/s41564-026-02496-z ; PRJNA943687 ; https://github.com/qsu123/ASD_multi-kingdom_diagnosis |
| C08 | Yap 2021, autism-related dietary preferences and microbiome associations; confounding context | https://doi.org/10.1016/j.cell.2021.10.015 ; use publisher's corrected version |
| C09 | Yano 2015, microbiota regulation of host serotonin biosynthesis | https://doi.org/10.1016/j.cell.2015.02.047 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC4393509/ |
| C10 | Williams 2014, characterized microbial tryptophan decarboxylases | https://pmc.ncbi.nlm.nih.gov/articles/PMC4260654/ |
| C11 | Host melatonin modulation and distinct tryptamine-derived receptor-active compounds | https://doi.org/10.1080/19490976.2024.2313769 ; https://doi.org/10.1021/acschembio.5c00313 |
| C12 | Zhao, bile-acid-associated IBS-D subset and experimental mechanisms | https://doi.org/10.1172/JCI130976 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC6934182/ |
| C13 | De Palma 2022, microbial histamine and strain-dependent visceral-pain mechanisms | https://doi.org/10.1126/scitranslmed.abj1895 ; https://pubmed.ncbi.nlm.nih.gov/35895832/ |
| C14 | 2024 fructan/gluten crossover and site-specific wheat-sensitivity work | https://doi.org/10.1186/s12916-024-03562-1 ; https://doi.org/10.1080/29933935.2024.2438621 |
| C15 | 2024 Hashimoto stool-metagenome/host-transcriptome study; 31 cases/30 controls, all female | https://doi.org/10.1186/s12967-024-05876-3 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC11577717/ |
| C16 | 2025 adult subclinical-hypothyroidism shotgun association study | https://doi.org/10.3390/microorganisms13112643 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC12654992/ |
| C17 | Feehley 2019, infant-microbiota transfer and food-allergy protection experiments | https://doi.org/10.1038/s41591-018-0324-z ; https://pmc.ncbi.nlm.nih.gov/articles/PMC6408964/ |
| C18 | Abdel-Gadir 2019, defined microbiota and immune-tolerance mechanisms | https://doi.org/10.1038/s41591-019-0461-z ; https://pmc.ncbi.nlm.nih.gov/articles/PMC6677395/ ; PRJNA525231 is 16S, not a shotgun training cohort |
| C19 | 2024 Anaerostipes caccae/lactulose synbiotic food-allergy experiments in mice | https://pmc.ncbi.nlm.nih.gov/articles/PMC11239278/ ; https://pubmed.ncbi.nlm.nih.gov/38906158/ |

### 15.4 Intervention, modelling and benchmark resources

| ID | Resource | Links and implementation availability |
|---|---|---|
| R01 | Deehan 2020, structure-specific resistant-starch intervention | https://doi.org/10.1016/j.chom.2020.01.006 ; https://pubmed.ncbi.nlm.nih.gov/32004499/ |
| R02 | Wastyk 2021, dietary fibre versus fermented foods | https://doi.org/10.1016/j.cell.2021.06.019 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC9020749/ |
| R03 | Lancaster 2022, distinct isolated dietary-fibre responses | https://pubmed.ncbi.nlm.nih.gov/35483363/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC9187607/ |
| R04 | Sardá 2025, baseline-dependent unripe-banana/inulin response | https://www.nature.com/articles/s41522-025-00817-4 |
| R05 | Song 2025, fibre-response prediction in prediabetes | https://www.nature.com/articles/s41467-025-66498-x ; https://pmc.ncbi.nlm.nih.gov/articles/PMC12749054/ ; raw projects CRA016259, CRA016260, CRA026307; clinical metadata require approved access and have noncommercial restrictions |
| R06 | Connors 2026, active-learning fibre/community experiments | https://www.nature.com/articles/s41589-026-02272-4 ; https://github.com/VenturelliLab/Connors-Thompson_et_al_2026 ; source/supplementary experimental tables accompany the article |
| R07 | Quinn-Bohmann 2024, community-model SCFA predictions | https://www.nature.com/articles/s41564-024-01728-4 ; https://github.com/Gibbons-Lab/scfa_predictions ; public validation projects PRJNA937304, PRJNA640404, PRJNA939256, PRJNA1033794; some human data require agreements |
| R08 | Quinn-Bohmann 2026, prebiotic/probiotic modelling across human trials | https://doi.org/10.1371/journal.pbio.3003638 ; https://github.com/Gibbons-Lab/2024_probiotic_engraftment ; https://zenodo.org/records/18037976 ; use final peer-reviewed methods, not superseded preprint performance claims; some source trial data are request-only |
| R09 | Fibre-tolerance counterevidence in specific contexts | https://doi.org/10.1053/j.gastro.2022.09.034 ; https://doi.org/10.1038/s41586-022-05380-y ; retain ex-vivo/animal model and population distinctions |
| D01 | Sylph primary methods and publicly described mocks | https://www.nature.com/articles/s41587-024-02412-y ; PRJEB52977 ; https://forgemia.inra.fr/metagenopolis/benchmark_mock |
| D02 | CAMI benchmarking resource | https://academic.oup.com/nar/article/53/W1/W102/8126258 |
| D03 | 2025 shotgun infectious-gastroenteritis validation | https://doi.org/10.1371/journal.pone.0331288 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC12404398/ |
| D04 | AMRFinderPlus primary methods and curation | https://pmc.ncbi.nlm.nih.gov/articles/PMC8208984/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC9455714/ |
| D05 | NCBI AMRFinderPlus implementation/database documentation | https://github.com/ncbi/amr |
| D06 | CARD/RGI curated resistance ontology and sequence models | https://card.mcmaster.ca/ ; pin database and tool releases independently |
| D07 | MEROPS peptidases, substrates and experimentally characterized specificity | https://www.ebi.ac.uk/merops/ ; use bacterial proteins and retain family/sequence provenance |

R08 is an additional validation resource for the optional model-assisted layer of the required synbiotic feature. It does not authorize changing the production donor matcher or claiming universal probiotic engraftment predictions. Reproduce its exact intervention/model context and distinguish publicly downloadable code from request-only trial inputs.

### 15.5 Data-ingestion checklist tied to features

| Asset family | Extract/store | Used by |
|---|---|---|
| dbCAN/CAZy/PUL | Protein/HMM assets, substrate mappings, experimentally characterized entries, CGC/PUL evidence | A01, mucin detail |
| Vitamin reconstructions | Required steps, valid alternatives, salvage/transport, positive/negative genomes | A03 |
| Reaction-specific sequence panels | Full sequences/accessions, diagnostic residues, negative homologs, locus context | A02/A04/A05/A06/A07/A08 |
| GBM definitions | Native mixed identifier namespaces and parser grammar, licensing, mapped reactions | Neuroactive/fermentation companion views |
| HydDB/BileActome | Sequence models, class definitions, training context, human-gut validation fixtures | Gas and bile transformation panels |
| Trait resources | Culture versus predicted evidence, strain/species mapping, uncertainty | Ecology/organism descriptions |
| Cohort resources | Participant/study IDs, eligibility, assay/tool versions, raw accession and outcome availability | New references and validation |
| Primary intervention records | Exact entity, target, model, endpoint, exposure, formulation, positive/null/conflicting findings | Action planner, synbiotic evidence matching |
| Strain–substrate records | SYN pairs, exact identities/preparations, measured growth versus predicted utilization, study arms and bundle outcomes | Required A16 evidence mode |
| Simulation resources | Model library, code version, solver, media/scenario constraints, public validation examples | A16 |

### 15.6 Personalized synbiotic research and downloadable assets

These records supplement R06–R08 and the intervention registry. All inputs needed for the new feature are identified here or in the embedded SYN table; no external handoff document is required. Retrieval and registry-building are implementation tasks. A reachable URL is not an installed, validated database. Pin the downloaded bytes and per-asset terms in the application manifest.

| ID | Source and direct assets | Required use and limits |
|---|---|---|
| S01 | ISAPP synbiotic consensus: https://www.nature.com/articles/s41575-020-0344-2 | Terminology and complementary/synergistic criteria; distinguish design intent from demonstrated host benefit. |
| S02 | Arzamasov 2025: https://doi.org/10.1038/s41564-025-02056-x ; raw growth: https://pmc.ncbi.nlm.nih.gov/articles/instance/12313528/bin/41564_2025_2056_MOESM7_ESM.xlsx ; metadata/substrates/validation: https://pmc.ncbi.nlm.nih.gov/articles/instance/12313528/bin/41564_2025_2056_MOESM4_ESM.xlsx ; references: https://doi.org/10.6084/m9.figshare.26053936 ; genomes PRJNA1126848 ; https://github.com/Arzamasov/compendium_manuscript ; https://github.com/Arzamasov/glycobif | Import measured growth separately from predicted pathways; store substrate identity and assay context. Pin code and archive versions during ingestion. |
| S03 | dbCAN-PUL: https://academic.oup.com/nar/article/49/D1/D523/5907960 ; https://pro.unl.edu/static/DBCAN-PUL/ ; https://pro.unl.edu/dbCAN_PUL/dbCAN_PUL/help | Public directory includes `dbCAN-PUL_Feb-2025.xlsx`, older versioned workbooks and sequence assets. Choose a compatible full snapshot; confirm terms for the exact files. Experimental locus evidence does not make every homolog a proven substrate user. |
| S04 | BacDive: https://academic.oup.com/nar/article/53/D1/D748/7848838 ; https://hub.dsmz.de/wiki/bacdive/api/ ; example https://api.bacdive.dsmz.de/v2/fetch/5621 ; terms https://bacdive.dsmz.de/about | API v2 is public without authentication as checked September 2026; retain measured/predicted distinction and collection aliases. Preserve CC-BY-4.0 and the site's additional commercial-use notice in terms metadata; do not flatten different notices into a software license. |
| S05 | Probio-Ichnos: https://www.mdpi.com/2076-2607/12/10/1955 ; https://github.com/Mtsif/probio-ichnos ; https://raw.githubusercontent.com/Mtsif/probio-ichnos/main/probioIchnos.json | Download checked: 3,265,140 bytes, 12,962 records, SHA256 `d352ac33d84e8dd127906eb27a430b12910c26483339e4eed2275c8ff01545a4`. Repo MIT. Contains PMID/species/strain and phenotype evidence; no genome-accession fields. Resolve publications and identity rather than treating scraped fields as validated clinical outcomes. |
| S06 | Integrated Probiotic Database: https://probiogenomics.unipr.it/index.php?route=/tools ; https://probiogenomics.unipr.it/files/IPDB_latest.zip ; https://www.oaepublish.com/articles/mrr.2024.11 | ZIP endpoint and header verified, announced size 100,125,954 bytes; full archive ingestion remains an implementation task. Mutable filename requires content hash and extracted inventory. Confirm archive terms independently of article terms. Optional genome discovery, not efficacy/availability evidence. |
| S07 | NCBI identity and downloads: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/genomes/download-genome/ ; https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/using-dataformat/genome-data-reports/ ; https://www.ncbi.nlm.nih.gov/biosample/docs/ | Resolve exact culture/BioSample/assembly aliases; retain versioned GCA/GCF and nucleotide records, not first species-name search hit. Download genomes/annotations/proteins with manifests and hashes. |
| S08 | Anaerostipes hadrus strain niches: https://pmc.ncbi.nlm.nih.gov/articles/PMC12503163/ ; growth https://pmc.ncbi.nlm.nih.gov/articles/instance/12503163/bin/table_s4_ycaf163.xlsx ; strain details https://pmc.ncbi.nlm.nih.gov/articles/instance/12503163/bin/table_s1_ycaf163.docx | Additional experimentally tested carbon-source edges; strain variation matters. Not clinical response evidence. |
| S09 | van Zanten 2012: https://doi.org/10.1371/journal.pone.0047212 ; growth figures https://pmc.ncbi.nlm.nih.gov/articles/instance/3474826/bin/pone.0047212.s001.tif | NCFM/Bl-04 carbohydrate screening and selected colon-culture pairs. Import tested outcomes, not merely a list of candidate carbohydrates. |
| S10 | R08 model code: https://github.com/Gibbons-Lab/2024_probiotic_engraftment/tree/a494f4e21c8510cf02618af1861f8649d92fe698 ; archive https://zenodo.org/records/18037976 | Archive metadata CC-BY-4.0; 21,628,338-byte ZIP MD5 `f29ee44aae961d76edcf6d3dfd955380`. Pin source commit and each data file. Source discrepancies and reproduction recipe are in section 11.7. |
| S11 | AGORA2 species model: https://zenodo.org/records/7739096/files/agora201_refseq216_species_1.qza?download=1 ; model method https://www.nature.com/articles/s41587-022-01628-0 ; workflows https://github.com/micom-dev/databases | Model archive CC-BY-SA-4.0, MD5 `6a13da2a3fd9b1bc059987d6344f32b6`, approximately 264.4 MB. Workflow Apache-2.0 differs from model-data terms. These species models do not establish exact commercial-strain behavior. |
| S12 | R06 culture design: https://doi.org/10.1038/s41589-026-02272-4 ; https://github.com/VenturelliLab/Connors-Thompson_et_al_2026/tree/b2fd741027a8243a59aad6bbd1883aed935e60d7 ; code https://doi.org/10.5281/zenodo.20336398 ; experimental data https://doi.org/10.5281/zenodo.19210709 | Archived code v1.0.1 and experiment metadata CC-BY-4.0. Synthetic/fecal projects PRJNA1233340/PRJNA1391539. Optional culture-domain reproduction; not a human personalized dosing engine. |
| S13 | Personalized-formulation precedent: https://journals.asm.org/doi/10.1128/msystems.00503-24 ; formulation supplement `msystems.00503-24-s0002.csv` on article | Open-label study with paired microbiome data and formulation records; allocation algorithm is proprietary. Do not label it an open trained recommender or causal validation of personalization. Published ingredients can inform curation within applicable reuse terms. |
| S14 | DS-01 whole-formulation trial: https://doi.org/10.3390/nu18020255 ; https://www.mdpi.com/2072-6643/18/2/255 ; NCT06009614 | Six-week adult placebo trial; 350 baseline, 219 evaluable at week six. Gastrointestinal outcomes support the tested bundle, not individual components or microbiome-guided allocation. Preserve attrition, funding and request-only data status. |
| S15 | Peng 2025 consortium/substrate screening: https://doi.org/10.1080/19490976.2025.2541028 ; supplement https://pmc.ncbi.nlm.nih.gov/articles/instance/12326573/bin/KGMI_A_2541028_SM0243.zip ; PRJNA1211324/PRJCA033703 | 73-strain consortium (44 commensals + 29 potential probiotics) with 28 edible substrates; not a 44×28 monoculture matrix. Some substrates are minerals, not established prebiotics. Human arm tested tablets alone. Preserve these distinctions and verify current raw-data access. |

**Exact S10 file downloads and SHA256 pins.** Prefix all four relative paths below with `https://raw.githubusercontent.com/Gibbons-Lab/2024_probiotic_engraftment/a494f4e21c8510cf02618af1861f8649d92fe698/`:

| Relative path | SHA256 of downloaded file | Use |
|---|---|---|
| `figures/source_data/S5_Data.xlsx` | `19ad4e8e2eea635f3c577df0d1bc04ba852249f368bac25d2519abbe71bd8c3e` | Released-output replay; paper supplement https://doi.org/10.1371/journal.pbio.3003638.s010 |
| `data/arivale-metagenomes.csv` | `e811818f16a1067a4e8d365a047bb0dbdf4c0fc3d185b89e840a53c3dd0edeea` | Public example input, 154 unique sample IDs |
| `european_medium.csv` | `0299fdd9b9d1ba2938681a6b078d54ddde73add1912a6b3582c8c9dffbea8f61` | Frozen baseline-medium input |
| `high_fiber_medium.csv` | `b29a95fc3bf1639b78bca65ea32d78ee7498753aa37ad653d366df9e2198d526` | Alternative-medium sensitivity |

Source entrypoints under the same prefix are `notebooks/arivale_probiotics.ipynb`, `notebooks/studyA_MCMMs.ipynb` and `notebooks/studyA_prediction.ipynb`. Inspect them with the explicit corrections in section 11.7. R08 study B raw data `PRJNA755324` and Arivale shotgun `PRJNA1262070` are additional validation inputs; study A shotgun inputs are request-only. A data-access request is not a required step for core evidence matching.

**S12 adapter details:** use `armored/models.py`, `armored/preprocessing.py` and `MiRNN Examples.ipynb` at the pinned commit. The model has 15 strain inputs and six fibre controls; outputs include 15 species-abundance trajectories plus four pH/metabolite trajectories: pH, lactate, butyrate and acetate. Training/solver dependencies and fitted weights are not supplied as a locked production bundle. Reproduce a small public culture example before a design search. Public strains and outcomes are downloadable at:

- https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41589-026-02272-4/MediaObjects/41589_2026_2272_MOESM4_ESM.xlsx
- https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41589-026-02272-4/MediaObjects/41589_2026_2272_MOESM10_ESM.xlsx

A useful experimental motif is *A. caccae* L1-92/DSM14662 + *B. uniformis* VPI0061/DSM6597 + inulin; *P. copri* CB7/DSM18205 modified its effect in that culture system. Do not substitute LAHUC for L1-92 or infer disease causation from this interaction. The human-derived validation uses fecal cultures, not human administration. Import the full strain table for broader candidate designs rather than replacing strain IDs with species names.

**Additional optional coverage:** APOLLO models, https://doi.org/10.7910/DVN/RL7E7G and https://www.sciencedirect.com/science/article/pii/S2405471225000298, may expand organism-model coverage, but exact bundle access/terms and interoperability were not verified here. Do not make this dataset mandatory. Retail label lookup may use https://dsld.od.nih.gov/api-guide ; a label is not independent verification of strain identity, viability or health benefit.

### 15.7 Research and data readiness across the entire addition

The research specifies how to implement the features, but **not every final accession registry, reference cohort or trained model is already assembled and validated**. This document is a complete work specification, not a claim of completed application validation. The implementation must close the tasks below and record actual state per asset and feature.

Use separate fields for `source_identified`, `retrievable`, `bytes_verified`, `registry_validated`, `example_reproduced` and `application_validated`. Do not collapse them into one “available” boolean. Keep `not_implemented`, `missing_input`, `access_restricted`, `scientifically_unresolved`, `not_applicable` and `computed` distinct. An engineering TODO cannot be disguised as a scientific limitation. Evidence cards do not satisfy a required quantitative feature unless that feature explicitly calls for evidence/context output.

| Capability | Available basis | Required closure before calling the new numerical feature validated |
|---|---|---|
| A01 substrates | dbCAN/PUL and primary substrate studies | Build all 14 exact substrate registries, physical/preparation specificity, positive/near-negative references, detection thresholds and compatible new calibration if percentiles are shown. |
| A02 fermentation | Explicit route logic, primary enzymes and HydDB | Resolve accession-level alternatives, direction, hydrogenase classes and competing reactions; quantify and test them. |
| A03 vitamins | Public synthesis/salvage tables and test genomes | Import a complete versioned reaction graph with alternative pathways and verified positive/negative genome fixtures. |
| A04 protein/nitrogen | MEROPS and reaction studies | Close the finite requested route inventory, reaction/sequence specificity and near-negative fixtures. |
| A05 neuroactive | GABA observations, GBM definitions, strain/consortium evidence | Resolve mixed historical IDs and sequence specificity. F17 now has an accession trail but not a validated serotonin-specific detector; implement the specified card and candidate-enzyme evidence without inventing that specificity. |
| A06 plant compounds | Exact ucdCFO anchor, glucosinolate locus and equol studies | Complete reaction-specific accession sets, regioselective alternatives and difficult negatives. |
| A07 GLP-1 mechanisms | P9 exact sequence/hash and primary evidence | Build/test detection against related proteases; host concentration remains a separately imported assay. |
| A08 mucin/LPS/bile | Public enzyme tables and BileActome | Import corrected sequence/substrate mappings, reaction alternatives and structural-evidence rules; validate on appropriate human data. |
| A09 ecology | Explicit descriptive formulas, trait resources and model methods | Freeze trait map and compatible cohort; train/freeze new community medoids/HACK adapter and any new reference distributions. Raw descriptive metrics can work before these model artifacts exist. |
| A10 explorer | Taxonomy fixtures and target inventories | Inspect installed application references, implement release-specific splits/aliases and verify claimed named-strain assessability. |
| A11 AMR | Standard databases and deduplication rules | Pin eligible database/tool versions; adapt actual current inputs and verify class/carrier/abundance semantics. |
| A12 actions | Existing library plus 60 I-records and the SYN records | Normalize exact preparation, endpoint, model, comparator and triggers; finish explicitly abstract-only extraction and deduplicate existing evidence. |
| A13 trends/context | Explicit formulas and comparability contract | Connect real specimen dates/hashes and historical objects; sufficient repeated measurements are input-dependent. |
| A14 inputs/imports | Complete schema/default/provenance contract | Inspect actual metadata graph; implement unit, identity, intake and lab adapters and consumers. |
| A15 layout | Placement/navigation contract | Integrate real templates; verify rendered report, links, numbering and unchanged canonical values. |
| A16 synbiotics | Exact pair seeds, public growth data, model code and public replay | Complete evidence-mode registry/ranker tests; for simulations install/pin compatible model/solver/media and reproduce a solved example. Personalized clinical response remains unvalidated. |

Cross-cutting work includes the section-3 accession registry and near-negative sets, new-method reference processing where required, asset-level terms/checksums, analytical validation on known positives/negatives/dilutions, actual repository integration and rendering. No new reference percentile is fabricated while its cohort is being built; show the real raw measurement and precise reference state.

For the ABS mechanism, the public `ABS_CHsu` repository at `aa7087436fb49f08d4e82c645ff3fc7640ab5827` lacks `bacteria_adhE_gene_labeled_seqs.fna`, `adhE_gene_genome_metadata.txt` and the required full cohort metadata. The publication's supplementary workbook does not supply the missing accession panel. Use an independently curated, explicitly named ethanol-route panel if the exact authors' assets remain unavailable; do not claim exact paper-panel reproduction. This is not permission to omit other ethanol-pathway capabilities.


## 16. Complete feature, taxonomy and action coverage manifest

This appendix is normative for completeness. IDs below identify report concepts and views, not independent clinical measurements. Bind each row to the actual protected legacy metric or to the named additive module. Repeated views must share a canonical underlying measurement while retaining their distinct context. Every new numerical output must use an explicit, versioned formula and documented reference definition. Existing scores remain unchanged.

### 16.1 Indexed measurement and interpretation concepts

| View ID | Concept | Binding / implementation |
|---|---|---|
| M001 | Protein breakdown | A04: proteolysis and protein/amino-acid fermentation capacities |
| M002 | Trimethylamine | Legacy TMA pathway; A15 precise TMA versus host TMAO labels |
| M003 | Ammonia | Legacy urease retained; A04 adds non-urease ammonia routes |
| M004 | Branched Chain Amino Acids | Legacy BCAA capacity; A04 overview |
| M005 | p-Cresol | Legacy p-cresol capacity |
| M006 | Histamine index | Legacy histamine capacity retained; A05/A15 expose components, not an undocumented replacement index |
| M007 | Methanobrevibacter smithii | Legacy methanogen detection; A10 archaeal taxon view |
| M008 | Methane production capacity | Legacy methane capacity |
| M009 | Beta-glucuronidase capacity | Legacy beta-glucuronidase capacity; A08/A14 distinguish DNA from measured activity |
| M010 | GABA breakdown | A05: separate GABA degradation; surface existing gabT evidence where compatible |
| M011 | Secondary bile acids | Legacy secondary-bile-acid result; A08 additional transformations |
| M012 | Antibiotic-resistance abundance index | A11: class-level determinant abundance; optional named ecological taxon surrogate remains distinct |
| M013 | Antibiotic-resistance richness index | A11: determinant richness; optional named ecological taxon surrogate remains distinct |
| M014 | Enterobacteriaceae | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M015 | Klebsiella | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M016 | Klebsiella pneumoniae | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M017 | Klebsiella oxytoca | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M018 | Salmonella enterica | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M019 | Escherichia coli | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M020 | Escherichia flexneri | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M021 | Escherichia dysenteriae | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M022 | Citrobacter | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M023 | Enterobacter | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M024 | Morganella | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M025 | Raoultella | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M026 | Streptococcus | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M027 | Staphylococcus | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M028 | Pseudomonas aeruginosa | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M029 | Haemophilus influenzae | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M030 | Haemophilus parainfluenzae | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M031 | Enterococcus faecium | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M032 | Enterococcus faecalis | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M033 | Clostridioides difficile | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M034 | Clostridium perfringens | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M035 | Acinetobacter baumannii | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M036 | Campylobacter | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M037 | Helicobacter pylori | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M038 | Blastocystis | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M039 | Cryptosporidium | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M040 | Entamoeba histolytica | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M041 | Entamoeba dispar | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M042 | Giardia | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M043 | Yersinia enterocolitica | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M044 | Vibrio | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M045 | Vibrio cholerae | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M046 | Cyclospora cayetanensis | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M047 | Candida | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M048 | Aspergillus | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M049 | Cryptococcus | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M050 | Saccharomyces | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M051 | Rhodotorula | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M052 | Saprochaete | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M053 | Malassezia | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M054 | Microsporum | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M055 | Trichophyton | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M056 | Bacteroidota | Legacy phylum composition; A10 hierarchical rollup |
| M057 | Firmicutes | Legacy phylum composition; A10 hierarchical rollup |
| M058 | Actinobacteriota | Legacy phylum composition; A10 alias-aware rollup |
| M059 | Proteobacteria | Legacy phylum composition; A10 alias-aware rollup |
| M060 | Bacteroides | A10 genus abundance with explicit member set |
| M061 | Bacteroides fragilis | Legacy organism result; A10 strain/toxin-qualified species view |
| M062 | Prevotella | A10 genus abundance with explicit member set |
| M063 | Ruminococcus | A10 historical-name-aware genus/group rollup |
| M064 | Ruminococcus gnavus | Legacy organism result; A10 descriptive view; shares measurement with M090 |
| M065 | Blautia | A10 genus/group abundance |
| M066 | Roseburia | A10 genus/group abundance |
| M067 | Phocaeicola dorei | Legacy organism result; A10 alias-aware view |
| M068 | Bifidobacterium | Legacy bifidobacterial measurements; A10 genus view; shares measurement with M073 |
| M069 | Faecalibacterium | Legacy organism group; A10 complete genus and species components |
| M070 | Akkermansia | Legacy organism result; A10 genus/species components |
| M071 | Firmicutes / Bacteroidota ratio | Legacy Firmicutes:Bacteroidota ratio unchanged; A09 companion dashboard |
| M072 | Proteobacteria / Actinobacteriota ratio | A09: descriptive Proteobacteria:Actinobacteriota ratio |
| M073 | Bifidobacterium | A10 probiotic-context genus view; same underlying total as M068 |
| M074 | Lactobacillaceae | A10 family rollup with explicit taxonomy release |
| M075 | Lacticaseibacillus rhamnosus | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M076 | Limosilactobacillus reuteri | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M077 | Lacticaseibacillus paracasei | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M078 | Lactiplantibacillus plantarum | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M079 | Lactobacillus acidophilus | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M080 | Bifidobacterium infantis | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M081 | Bifidobacterium bifidum | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M082 | Bifidobacterium longum | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M083 | Bifidobacterium breve | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M084 | Bifidobacterium animalis | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M085 | Bifidobacterium adolescentis | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M086 | Streptococcus thermophilus | Legacy detection retained; A10 explicit taxon screening/rollup view with aliases, assay coverage and strain context |
| M087 | Shannon diversity | Legacy Shannon preserved; A09 states catalogue, filtering and logarithm base |
| M088 | Species richness | Legacy total and reference-comparable richness preserved; A09/A10 show both |
| M089 | Overabundant species | Legacy organism expansion retained; A09 top-organism dominance view |
| M090 | Ruminococcus gnavus | Legacy expansion interpretation; A10 shares organism measurement with M064 |
| M091 | Gut resilience score | A09: experimental ecological resilience potential plus measured longitudinal recovery |
| M092 | Oral microbes | Legacy oral-origin group; A10 components and exclusions |
| M094 | Butyrate | Legacy butyrate capacity; A02 network context |
| M095 | Propionate | Legacy propionate capacity; A02 network context |
| M096 | Acetate | A02: acetate routes |
| M097 | Cellulose | A01: cellulose utilization |
| M098 | Resistant starch | A01: resistant-starch utilization |
| M099 | Chitin | A01: chitin utilization |
| M100 | Pectin | A01: pectin utilization |
| M101 | Fructooligosaccharides | A01: fructooligosaccharide utilization |
| M102 | Galactooligosaccharides | A01: galactooligosaccharide utilization |
| M103 | Xylooligosaccharides | A01: xylooligosaccharide utilization |
| M104 | Isomaltooligosaccharides | A01: isomaltooligosaccharide utilization |
| M105 | Vitamin B2 | Legacy riboflavin capacity; A03 complete vitamin overview |
| M106 | Vitamin B7 | Legacy biotin capacity; A03 overview |
| M107 | Vitamin B9 | Legacy folate capacity; A03 overview |
| M108 | Vitamin B12 | Legacy cobalamin capacity; A03 synthesis/salvage distinction |
| M109 | Vitamin K | Legacy menaquinone capacity; A03 distinguishes vitamin forms |
| M110 | Indole-3-propionic acid | Legacy indole-3-propionate capacity |
| M111 | GABA production | Protected legacy GABA production; A05 adds separate degradation view |
| M112 | Unconjugated bile acids | Legacy bile-salt hydrolase/unconjugated-bile-acid result; A08 detail |
| M113 | Urolithin-producing species | Legacy polyphenol guild retained; A06 separate urolithin contributors and reactions |
| M114 | Oxalate degradation | Legacy oxalate capacity |
| M115 | Hexa-LPS index | Legacy hexa-acylated-LPS carrier result; A08 additional structural-gene evidence |
| M116 | Mucus degradation index | Legacy mucin-degrader group retained; A08 substrate-specific enzymes |
| M117 | Hydrogen sulfide index | Legacy H2S capacity; A02/A04 source-route context |
| M118 | Host DNA | A14/A13 provider host-DNA QC import when original prefilter denominator exists |
| M119 | Oxygen exposure index | A09: ecological aerotolerance proxy; no claim of directly measured oxygen leakage |

**Explicit repeats:** M064/M090 are the same organism in descriptive and expansion contexts; M068/M073 are the same genus total in resident/beneficial and probiotic contexts. Do not count them as four independent assays. Genus/family totals must document any excluded taxa; a label implying a complete genus cannot silently mean an opportunist-only subset.

### 16.2 Extended functional-concept coverage

Every functional concept below must resolve, including those already represented above. A mechanism/evidence card or measured external assay is the correct implementation where microbial DNA does not identify a host metabolite concentration.

| Function ID | Concept | Binding / implementation |
|---|---|---|
| F001 | Acetate | M096 / A02 |
| F002 | Butyrate | M094 / protected legacy; A02 context |
| F003 | Propionate | M095 / protected legacy; A02 context |
| F004 | Acetylcholine | A05 evidence and route-specific microbial capacity where demonstrated; A14 measured endpoint import |
| F005 | Dopamine | Protected legacy tyrosine-decarboxylase result; A05 reaction-specific catecholamine handling, not brain dopamine concentration |
| F006 | GABA | M111 production plus M010 degradation; A05 |
| F007 | Histamine | M006 / protected legacy; A05 |
| F008 | Norepinephrine | A05 mechanistic evidence/validated reaction only; A14 measured endpoint import |
| F009 | Serotonin | A05 experimentally supported microbial routes and host-signalling context; A14 measured endpoint import |
| F010 | Vitamin B1 (Thiamin) | A03 thiamine synthesis/salvage |
| F011 | Vitamin B2 (Riboflavin) | M105 / protected legacy; A03 |
| F012 | Vitamin B3 (Niacin) | A03 niacin/NAD precursor synthesis/salvage |
| F013 | Vitamin B5 (Pantothenic Acid) | A03 pantothenate synthesis/salvage |
| F014 | Vitamin B6 (Pyridoxine) | A03 B6-vitamer synthesis/salvage |
| F015 | Vitamin B7 (Biotin) | M106 / protected legacy; A03 |
| F016 | Vitamin B9 (Folate) | M107 / protected legacy; A03 |
| F017 | Vitamin B12 (Cobalamin) | M108 / protected legacy; A03 |
| F018 | Vitamin K2 | M109 / protected legacy; A03 |
| F019 | Ammonia | M003 / legacy urease plus A04 additional routes |
| F020 | Hydrogen sulfide | M117 / protected legacy; source-specific routes |
| F021 | Lipopolysaccharides (LPS) | M115 / legacy carrier view plus A08 lipid-A structural potential |
| F022 | Methane | M008 / protected legacy |
| F023 | Trimethylamine (TMA) | M002 / protected legacy; preserve host TMAO distinction |
| F024 | Beta-glucuronidase | M009 / legacy beta-glucuronidase capacity; A14 measured activity |
| F025 | GLP-1 activators | A07 component-based microbial mechanism panel; A14 measured host GLP-1 |
| F026 | Indoles and Phenols | Preserve separate legacy indole, IPA and p-cresol results; A04/A05 grouped navigation without conflating products |
| F027 | Myrosinase | A06 glucosinolate-to-isothiocyanate conversion mechanisms, distinguishing microbial routes from plant enzyme |
| F028 | Urolithin | M113 / A06 exact conversion reactions and taxa |
| F029 | Chitin | M099 / A01 |
| F030 | Resistant starch | M098 / A01 |
| F031 | Pectin | M100 / A01 |
| F032 | Inulin | A01 fructan utilization with degree-of-polymerization context |
| F033 | Fructooligosaccharides | M101 / A01 |
| F034 | Xylooligosaccharides | M103 / A01 |
| F035 | Galactooligosaccharides | M102 / A01 |
| F036 | Protein metabolizers | M001 / A04; taxon guild distinguished from functional genes |
| F037 | Lactose degraders | A01 microbial lactose utilization; does not measure human lactase deficiency |
| F038 | Oxalate degraders | M114 / protected legacy; distinguish taxa from pathway support |

### 16.3 Context and symptom navigation coverage

A13 must expose links from supported findings to these contexts: anxiety/depression; autism research; serotonin/melatonin-related gut physiology; atherosclerosis; hypertension; auto-brewery/ethanol-production evaluation; chronic constipation; chronic diarrhea; recurring abdominal pain; gluten-related symptoms; IBS; thyroid-related context; food allergy; GLP-1-related microbial mechanisms; eczema. Existing IBS and other disease-profile calculations already present are preserved, not reimplemented as new diagnoses. Additional contexts may have an evidence card without a fabricated classifier. Every context lists all contributing features rather than hiding them behind a truncated count.

### 16.4 Taxonomy alias-resolution fixture

The **taxonomy labels below** form the taxonomy-resolution regression fixture for A10. They are not required positive detections in any specimen, a reference abundance vector, or assertions of health effects. Resolve each to a pinned taxonomy concept/cluster where evidence allows; otherwise retain the original string with an explicit unresolved or one-to-many mapping. Do not invent a one-to-one synonym, merge genuinely distinct clusters, or count aliases twice. All suffixes and provisional names are intentional.

```text
Phocaeicola vulgatus, Ruminococcus_B gnavus, Blautia_A massiliensis, Bacteroides uniformis,
Phocaeicola plebeius_A, Parabacteroides distasonis, Bacteroides stercoris, Anaerostipes sp900756035,
Alistipes putredinis, Fusicatenibacter saccharivorans, Bacteroides thetaiotaomicron, Bifidobacterium longum,
Bacteroides caccae, Alistipes onderdonkii, Barnesiella intestinihominis, Parabacteroides merdae,
Bifidobacterium breve, Bacteroides ovatus, Lactiplantibacillus plantarum, Paraprevotella clara,
Blautia_A caecimuris, Collinsella sp900545605, Phocaeicola sp900554435, UBA9502 sp900538475,
Phocaeicola massiliensis, Anaerotignum lactatifermentans, Akkermansia muciniphila, Enterocloster sp001517625,
Blautia sp000432195, Gemmiger sp900545545, Alistipes shahii, UBA9475 sp002161675,
Flavonifractor plautii, Blautia_A sp003471165, Gemmiger variabilis, Faecalimonas phoceensis,
Bifidobacterium infantis, Mediterraneibacter faecis, Odoribacter splanchnicus, Bifidobacterium animalis,
Eubacterium_I ramulus, Bariatricus comes, Bacteroides sp900066265, Dialister sp000434475,
Sutterella wadsworthensis, Duodenibacillus sp900538905, CAG-877 sp900548695, Erysipelatoclostridium ramosum,
UBA1691 sp900544375, Blautia_A wexlerae, Bacillus velezensis, Blautia_A sp900541345,
Bifidobacterium adolescentis, Schaedlerella glycyrrhizinilytica, Dorea formicigenerans, CAG-317 sp900543415,
Collinsella sp003438495, Ruthenibacterium lactatiformans, Enterocloster bolteae, Butyricimonas faecihominis,
Clostridium_Q saccharolyticum_A, Bacteroides xylanisolvens, Flavonifractor sp000508885, Eisenbergiella massiliensis,
Faecalibacterium sp900540455, Phocaeicola dorei, Gemmiger sp900540595, Blautia_A obeum,
Bacteroides sp014385165, Eisenbergiella sp900544445, Blautia sp003287895, Lacticaseibacillus rhamnosus,
Phocaeicola plebeius, Erysipelatoclostridium spiroforme, Lachnospira sp000437735, Fournierella massiliensis,
Ruminococcus_B sp900544395, Collinsella aerofaciens_G, Clostridium_A leptum, Evtepia gabavorous,
Coprobacter fastidiosus, Dysosmobacter welbionis, Collinsella sp900541235, Eggerthella lenta,
Bilophila wadsworthia, Senegalimassilia anaerobia, Lachnoclostridium_B sp900066555, Sellimonas intestinalis,
CAG-81 sp900066535, Collinsella aerofaciens_J, Streptococcus thermophilus, Lawsonibacter sp902363045,
Agathobacter rectalis, Gemmiger formicilis, Gemmiger sp900540775, Muricomes oroticus,
Ruminococcus_A sp000432335, Limosilactobacillus reuteri, Gordonibacter pamelaeae, Eubacterium maltosivorans,
CAG-317 sp000433535, Clostridium_AQ innocuum, Agathobaculum butyriciproducens, Faecalimonas umbilicata,
Lawsonibacter asaccharolyticus, CAG-273 sp000438355, Agathobaculum sp900291975, Ligilactobacillus salivarius,
Lawsonibacter sp900066825, Weizmannia coagulans_A, Blautia_A wexlerae_A, UMGS1670 sp900546215,
Enterocloster citroniae, Holdemania massiliensis, Oxalobacter sp900760095, Anaerofustis stercorihominis,
Blautia_A sp003480185, CAG-81 sp000435795, Mediterraneibacter torques, Lawsonibacter sp900066645,
Gemmiger variabilis_C, Acutalibacter sp900548545, Blautia sp900753905, Blautia sp002161285,
Borkfalkia ceftriaxoniphila, Lactobacillus paragasseri, Blautia_A sp900066205, Enterococcus_B faecium,
Anaerostipes hadrus_B, Alistipes communis, Enterocloster clostridioformis, OF09-33XD sp003481995,
Phocaeicola mediterraneensis, Anaerotignum sp001304995, Butyricimonas paravirosa, Enterocloster aldenensis,
Faecalibacterium prausnitzii_C, Negativibacillus massiliensis, CAG-81 sp900541255, Dysosmobacter sp014297375,
Anaeromassilibacillus sp002159845, Anaerofilum sp002160015, CAG-1427 sp000436075, Collinsella aerofaciens_F,
Dialister sp900543165, Butyricicoccus pullicaecorum, UBA9475 sp900549885, Anaerosacchariphilus sp900066385,
Collinsella sp900540985, Lachnospira sp900552795, Bacteroides sp900761785, UMGS1071 sp900542375,
An92 sp900199495, Delmidovirus intestinihominis, Duodenibacillus sp003472385, Hungatella effluvii,
Bacteroides finegoldii, Blautia sp001304935, Collinsella sp900751755, Dysosmobacter sp900548505,
Gemmiger qucibialis, Sellimonas sp002161525, Alistipes_A indistinctus, Collinsella sp900550205,
Lactobacillus acidophilus, UBA3402 sp003478355, Dialister sp900543455, Intestinimonas butyriciproducens,
Blautia hansenii, Collinsella sp900758475, Clostridium_Q sp003024715, Lawsonibacter sp000177015,
Parabacteroides sp900541965, Anaerotruncus colihominis, Blautia sp900547685, Enterocloster lavalensis,
Gordonibacter urolithinfaciens, Paraprevotella xylaniphila, Afonbuvirus faecalis, Blautia_A sp900066335,
Blautia_A sp900549015, Enterocloster sp900547035, Collinsella sp900754435, Dorea_A longicatena,
Anaerobutyricum sp900554965, Lacticaseibacillus paracasei, Gemmiger sp900539695, Gemmiger sp900554145,
Intestinimonas timonensis, Schaedlerella sp900066545, UBA1417 sp900552925, UBA9502 sp900540335,
BX12 sp014333425, CAG-83 sp000435975, Collinsella aerofaciens, Eisenbergiella sp900539715,
Faecalibacterium prausnitzii_G, Faecalimonas sp900550975, Marseille-P3106 sp900169975, UMGS775 sp900545325,
Bacteroides sp900765785, Blautia_A sp003474435, Acetatifactor sp900066565, Agathobaculum sp900557315,
Collinsella sp900758375, Hydrogeniiclostridium mannosilyticum, Phocaeicola sp900557085, Senegalimassilia sp900550055,
UBA1691 sp900544715, UBA6398 sp003150315, Anaerosacchariphilus sp002160825, Bacteroides cellulosilyticus,
Bifidobacterium bifidum, Blautia_A schinkii, Blautia_A sp900548245, Dorea_B phocaeensis,
Gemmiger_A sp002160955, HGM12998 sp900756495, AM51-8 sp900761925, Bacteroides sp003463205,
Collinsella sp900765115, Delmidovirus splanchnicus, Intestinimonas massiliensis, Blautia sp900120295,
Collinsella sp900541885, Collinsella sp900752015, Eisenbergiella sp900555195, Faecalibacterium prausnitzii,
Lawsonibacter sp900545895, Lawsonibacter sp900754605, Parabacteroides acidifaciens, Phocaeicola sartorii,
Streptococcus anginosus, UBA1191 sp900549125, UMGS966 sp900547185, CAG-145 sp000435715,
Collinsella sp003487125, Gemmiger variabilis_B, Lachnoclostridium_B sp900543315, Bacteroides fragilis_A,
Blautia sp900539145, CAG-145 sp900542565, CAG-83 sp900549395, Collinsella sp900542085,
Firm-11 sp900540045, Lawsonibacter sp002160305, Longicatena caecimuris, Mediterraneibacter massiliensis,
Pediococcus pentosaceus, Pseudoflavonifractor capillosus, QALS01 sp900552625, Ruthenibacterium sp003149955,
Bacteroides sp003545565, Bacteroides sp003865075, Clostridium_AQ sp003481775, Clostridium_Q sp000435655,
Collinsella sp002232035, Enterocloster sp000431375, Phocaeicola sp000436795, Phocaeicola sp002493165,
Agathobacter faecis, Angelakisella massiliensis, Bacteroides faecis, Bacteroides fragilis,
Blautia_A faecis, Blautia_A sp900553515, CAG-41 sp900066215, Clostridium_Q symbiosum,
Collinsella sp900544845, Collinsella sp900762015, Faecalibacterium prausnitzii_D, Faecalibacterium sp002160895,
Gemmiger variabilis_A, Holdemanella sp002299315, Massilimaliae massiliensis, UBA9475 sp002161235,
UMGS1537 sp900552695, Bacteroides eggerthii, Blautia_A sp900542045, Blautia_A sp900547615,
Butyricimonas virosa, CAG-81 sp900066785, Clostridium_AP scindens, Fournierella sp002161595,
Intestinibacter sp900553485, Lachnoclostridium_A edouardi, Lachnoclostridium_B phocaeensis, NK3B98 sp900758315,
Parasutterella excrementihominis, Alistipes finegoldii, Bittarella massiliensis, Collinsella sp003470665,
Collinsella sp900541145, Collinsella sp900759335, Dorea_A longicatena_B, Dysosmobacter sp001916835,
Eubacterium_G sp900556905, Frisingicoccus caecimuris, Holdemania filiformis, Lactiplantibacillus paraplantarum,
Lactonifactor sp009677585, Phocaeicola coprophilus, Prevotella copri, UBA1191 sp900545775,
Anaerobutyricum hallii, Anaerobutyricum soehngenii, Anaeromassilibacillus sp001305115, Anaerostipes hadrus,
Blautia sp900541955, Blautia_A sp900120195, Collinsella aerofaciens_M, Collinsella sp003436275,
Collinsella sp003439125, Collinsella sp900544645, Escherichia flexneri, Faecalibacterium prausnitzii_I,
Finegoldia magna, Finegoldia magna_H, Fournierella sp002159185, Fournierella sp004558145,
Lachnoclostridium_A sp002160755, Limosilactobacillus reuteri_E, Massilistercora timonensis, Muricomes contortus_B,
Muricomes sp900604355, Phocaeicola coprocola, Phocaeicola sp002161565, Propionibacterium freudenreichii,
Pseudoruminococcus massiliensis, Schaedlerella sp004556565, Thauera sp009469595, UBA1096 sp002684915
```

### 16.5 Canonical food and ingredient vocabulary

A12 must support the following vocabulary, including links from a broad food class to its specific members. These are evidence-retrieval and shopping-list terms, not automatic recommendations. Retain food form, processing, dose context and formulation separately; edible garlic, purified allicin, an essential oil and an extract are different interventions. Canonical IDs below use the `action.` prefix.

| Registry ID | Canonical concept and explicitly supported aliases/forms |
|---|---|
| `action.food.apple` | Apples |
| `action.food.asparagus` | Asparagus |
| `action.food.barley` | Barley |
| `action.food.olive` | Olives: black and green |
| `action.food.berry` | Berries; blueberry, blackberry, strawberry, raspberry as separate child terms |
| `action.food.cabbage` | Cabbage |
| `action.food.cacao` | Cacao/cocoa powder; cocoa husks remain a distinct preparation |
| `action.food.caper` | Capers |
| `action.spice.cardamom` | Cardamom |
| `action.food.cashew` | Cashew |
| `action.food.celery` | Celery; celery seed as a distinct spice preparation |
| `action.spice.cinnamon.ceylon` | Ceylon cinnamon |
| `action.food.chicory` | Chicory greens; chicory root; inulin extract as a distinct ingredient |
| `action.spice.clove` | Cloves |
| `action.food.cranberry` | Cranberries |
| `action.spice.cumin` | Cumin |
| `action.food.dairy` | Dairy; keep fermented/unfermented and product composition distinct |
| `action.food.dandelion` | Dandelion greens |
| `action.spice.mixture` | Spice mixtures with individually identified components; includes star anise and peppermint |
| `action.food.flax` | Flaxseed; flaxmeal |
| `action.spice.oregano` | Fresh oregano; extracts/oils are separate formulations |
| `action.food.fruit` | Fruit category; does not replace named fruit terms |
| `action.food.ginger` | Ginger; ginger juice as distinct preparation |
| `action.food.artichoke.globe` | Globe artichoke/green artichoke heads |
| `action.food.green_bean` | Green beans |
| `action.food.nut` | Nuts: hazelnut, pecan, almond, walnut, cashew, pistachio as child terms |
| `action.food.leek` | Leeks |
| `action.food.legume` | Legumes: beans, lentils and chickpeas; navy bean as a child term |
| `action.food.mushroom` | Mushrooms: shiitake, maitake and reishi; extracts separate from whole food |
| `action.food.yeast` | Nutritional yeast and brewer’s yeast; live versus inactive preparation |
| `action.food.seed` | Seeds: chia and flax; separate from nut allergy category |
| `action.food.oat` | Oats: rolled, steel-cut and oat bran as preparation-specific children |
| `action.food.onion` | Onions |
| `action.food.pear` | Pears |
| `action.spice.peppermint` | Peppermint leaf; oil/extract separate |
| `action.food.plum` | Plums |
| `action.food.pomegranate` | Pomegranate; extract separate |
| `action.food.garlic` | Raw garlic; cooked garlic, juice and purified compounds distinct |
| `action.food.rye` | Rye |
| `action.spice.sage` | Sage |
| `action.food.seaweed` | Seaweed: wakame and kelp, including Laminaria digitata; species and extract composition retained |
| `action.food.jerusalem_artichoke` | Sunchoke/Jerusalem artichoke; distinct from globe artichoke |
| `action.food.sweet_potato` | Sweet potato; resistant-starch exposure depends on preparation |
| `action.spice.thyme` | Thyme; thymol and thyme essential oil are separate interventions |
| `action.food.vegetable` | Vegetables category; retain named vegetables |
| `action.food.wheat` | Wheat; wheat bran and wheat flour as distinct preparations |
| `action.food.whole_grain` | Whole grains category; retain grain and preparation |
| `action.food.yacon` | Yacon syrup; composition/formulation retained |
| `action.food.annatto` | Annatto; identify preparation and compounds |
| `action.food.broccoli_sprout` | Broccoli sprouts; separate glucoraphanin, myrosinase and sulforaphane preparations |
| `action.food.rice` | Rice; retain grain type, processing and cooling preparation |
| `action.compound.gallic_acid` | Gallic acid; food additive or defined compound |
| `action.food.vegetable_oil_blend` | Vegetable oil blend; exact component oils and proportions required |
| `action.food.soy` | Soy; soybeans and soymilk as separate preparations |
| `action.food.chicken` | Chicken; preparation and comparator diet recorded |
| `action.food.kiwifruit` | Kiwifruit; cultivar/formulation retained |
| `action.food.konjac` | Konjac/konjaku flour; glucomannan as a defined ingredient |
| `action.food.rhubarb` | Rhubarb; species, plant part and preparation required |
| `action.food.agave_syrup` | Agave syrup; distinct from agave fructans/inulin |
| `action.vitamin.a` | Vitamin A; retinoid form and existing vitamin intake recorded |
| `action.compound.xylitol` | Xylitol; oral/dental and stool evidence kept separate |
| `action.mineral.iron` | Iron; salt, route and measured indication retained |
| `action.polyphenol.lychee` | Lychee/litchi-derived polyphenol preparation |
| `action.polyphenol.isoflavone` | Isoflavones; individual compounds and mixture composition retained |
| `action.compound.propolis` | Propolis/bee glue; composition and allergy context retained |
| `action.prebiotic.xos` | Xylooligosaccharides/XOS; degree of polymerization recorded |
| `action.prebiotic.beta_glucan` | Beta-glucans; cereal versus yeast/fungal linkage structure retained |
| `action.prebiotic.phgg` | Partially hydrolysed guar gum/PHGG; distinct from native guar |
| `action.prebiotic.arabinoxylan` | Arabinoxylan; source and hydrolysis retained |
| `action.prebiotic.inulin` | Inulin; chain length and formulation retained |
| `action.polyphenol.genistein` | Genistein; distinct from an unspecified isoflavone mixture |
| `action.polyphenol.tannin` | Tannin extracts; exact chemical mixture and plant source required |
| `action.botanical.pink_trumpet_tree` | Pink trumpet tree/Tabebuia impetiginosa preparation; botanical identity, plant part and extract required |
| `action.prebiotic.laminarin` | Laminarin-enriched seaweed extract; distinct from whole seaweed |
| `action.probiotic.lacticaseibacillus_rhamnosus_hn001` | Lacticaseibacillus rhamnosus HN001; historical Lactobacillus rhamnosus HN001 |
| `action.probiotic.enterococcus_faecium` | Enterococcus faecium preparation: species-only evidence remains unresolved until strain, resistance and product are identified |
| `action.probiotic.bifidobacterium_longum_bb536` | Bifidobacterium longum BB536 |
| `action.probiotic.bifidobacterium_animalis_lactis` | Bifidobacterium animalis subsp. lactis: strain required for strain-specific claims |
| `action.probiotic.enterococcus_durans_ep1` | Enterococcus durans EP1 |
| `action.postbiotic.eps_bifidobacterium_longum` | Purified exopolysaccharide from specified Bifidobacterium longum strain; not interchangeable with live species |
| `action.postbiotic.eps_bifidobacterium_animalis` | Purified exopolysaccharide from specified Bifidobacterium animalis strain |
| `action.postbiotic.eps_bifidobacterium_pseudocatenulatum` | Purified exopolysaccharide from specified Bifidobacterium pseudocatenulatum strain |
| `action.pattern.vegetarian` | Vegetarian/plant-based dietary pattern; composition retained |
| `action.pattern.gluten_free` | Gluten-free dietary pattern; diagnostic context and substitutions retained |
| `action.pattern.high_fibre` | High-fibre dietary pattern; substrate mixture and tolerance retained |
| `action.pattern.mediterranean` | Mediterranean dietary pattern |
| `action.lifestyle.exercise` | Exercise/physical activity; type, intensity and duration retained |
| `action.lifestyle.weight_change` | Weight-loss intervention; distinguish intervention from the resulting weight change |
| `action.pattern.ketogenic` | Ketogenic dietary pattern; composition and comparator retained |
| `action.lifestyle.cold_exposure` | Cold exposure; protocol and population retained; no assumed gut benefit |
| `action.supplement.butyrate` | Butyrate/tributyrin products with formulation and release-site metadata |
| `action.probiotic.bacillus` | Bacillus/related-genus probiotic products with exact species, strain and formulation |
| `action.product.synbiotic` | Synbiotic products with exact probiotic strains and prebiotic components |
| `action.supplement.polyphenol_mixture` | Polyphenol mixtures with defined ingredient composition |
| `action.prebiotic.fos` | Fructooligosaccharides/FOS |
| `action.prebiotic.cellulose` | Cellulose; preparation and co-administered food matrix |

Implement every canonical food and intervention concept in the vocabulary above; broad categories such as nuts/seeds, legumes, fruit and vegetables retain their named children. Merge duplicate shopping-list entries, not biologically different preparations. The existing action catalogue remains intact and can contain additional ingredients.

#### Target retrieval terms and evidence binding

Use these canonical organism concepts when searching and binding action evidence: `Akkermansia muciniphila`, `Bifidobacterium` genus, `Bifidobacterium longum`, `Blautia wexlerae` (plus supported current-taxonomy aliases), `Escherichia coli`, `Faecalibacterium` genus, `Faecalibacterium prausnitzii` complex, `Fusobacterium` genus, `Roseburia` genus, `Streptococcus thermophilus`, and `Enterococcus faecium`. Functional targets include acetate/SCFA potential, butyrate potential, TMA formation, beta-glucuronidase capacity, cellulose utilization, and resilience context. These are retrieval keys, not a request to increase or suppress every listed organism.

Each intervention-to-target edge requires an independently verified source, model/body site, exact preparation or strain, measured direction, endpoint, population, competing effects and confidence. An observation that a food reduces one taxon does not by itself make that food undesirable or establish a clinical benefit. Species-only probiotic evidence cannot authorize a strain-specific product claim. Preserve laboratory and animal findings as explicitly labeled experimental evidence; do not falsely promote them to human outcome evidence.

### 16.6 Body-site and assay exclusions

- Vaginal community-state typing, vaginal Lactobacillus dominance and vaginal infection panels require a vaginal specimen. They are not derived from stool.
- Oral/dental microbiome outcomes, including oral xylitol or Streptococcus effects, do not establish the same change in stool. Oral-origin organisms detected in stool remain a separate ecological finding.
- Skin/scalp fungal studies do not establish gut fungal overgrowth or treatment response. Species associated with skin or environmental exposure may still appear in stool screening with appropriate specimen context.
- Host hormones, inflammatory proteins, pancreatic function, measured stool pH, fecal blood, metabolite concentrations and enzymatic activity require the relevant external assay. Genetic capacity cards may link to them but must not supply invented values.
- Gastrointestinal sub-sites and stool are distinct. A stool H. pylori DNA result does not directly measure gastric mucosal inflammation; a stool taxon does not prove invasive infection elsewhere.
- DNA-only libraries do not generally establish RNA-virus absence. A target outside the assay’s molecular scope receives `unsupported_by_assay`, not a reassuring negative.

### 16.7 Completeness acceptance

The implementation manifest must contain every M-series row listed above and F001–F038, with explicit aliasing for repeated concepts, plus every taxonomy label in section 16.4 and the canonical intervention vocabulary. Render every implemented result in the complete index; every unavailable result must state its actual data requirement. Require complete implementation coverage of this embedded manifest. Validate numerical outputs against declared formulas, reference fixtures and independently established ground truth; validate organism calls against sequence evidence and analytical controls.

## 17. Intervention evidence seed matrix

All 60 source URLs resolved to the corresponding publication records. Titles, DOI identifiers, publication type and complete abstracts were retrieved from Europe PMC/MEDLINE on 22 September 2026 and manually reviewed. This is abstract-level verification; full text is still required before promoting a specific effect size, exact assay concentration, product strain or species direction not explicitly documented in the abstract. Two cited works are reviews, and one uses a case narrative; do not promote these as primary trials.

Study protocols below are source facts, not dosing recommendations. Human, animal, cell and ex-vivo evidence are deliberately included. A null or contrary result remains visible. Add these as candidate evidence cards in new v0.8.3 modules, preserving existing functionality.

Each I001–I060 ID identifies a source record. Derive separate intervention–target edges for distinct findings, using the schema in section 8; do not treat one table row as a universal prescription.

| Evidence ID | Exact ingredient / strain | Model | Supported endpoint | Applicability / limitation | Primary citation |
|---|---|---|---|---|---|
| I001 | Oats, 80 g/day versus rice, 45 days | Human randomized, 210 adults | Increased Akkermansia and Roseburia; improved lipid endpoints | Strong substrate/action seed; actual oats intervention, not all isolated oat fibers | [PMID 34956218](https://pubmed.ncbi.nlm.nih.gov/34956218/) · DOI 10.3389/fimmu.2021.787797 |
| I002 | Cooked navy bean powder, 15.7% diet | Mouse obesity experiment | Increased Akkermansia, Prevotella and SCFA; improved intestinal and adipose markers | Animal model; no human dose inference | [PMID 33652785](https://pubmed.ncbi.nlm.nih.gov/33652785/) · DOI 10.3390/nu13030757 |
| I003 | Annatto-extracted tocotrienol, 800 mg/kg diet | Mouse high-fat-diet experiment | Glucose and adipokine improvement; Akkermansia higher versus low-fat comparator | Extract differs from annatto food coloring; comparator-specific direction | [PMID 32438021](https://pubmed.ncbi.nlm.nih.gov/32438021/) · DOI 10.1016/j.nutres.2020.04.001 |
| I004 | Walnuts, 42 g/day | Human randomized crossover, 18 adults | Increased Faecalibacterium/Roseburia; lower fecal DCA/LCA; Bifidobacterium also lower | Show multiple effects and measured bile acid outcome; lowering bifidobacteria is not itself the therapeutic objective | [PMID 29726951](https://pubmed.ncbi.nlm.nih.gov/29726951/) · DOI 10.1093/jn/nxy004 |
| I005 | Raw broccoli sprouts, 20 g/day versus alfalfa sprouts | Human intervention, 48 adults with constipation scores | Improved bowel score and defecation duration; Bifidobacterium percentage decreased | Relevant glucosinolate/sulforaphane card; not proof that low microbial myrosinase warrants supplementation | [PMID 29371757](https://pubmed.ncbi.nlm.nih.gov/29371757/) · DOI 10.3164/jcbn.17-42 |
| I006 | Whole black raspberries, anthocyanin or residue fractions | Rat feeding experiment | Whole berries increased Akkermansia and Anaerostipes; fractions differed | Use black raspberry product/model; abstract does not support generic suppression of bifidobacteria as target | [PMID 28718724](https://pubmed.ncbi.nlm.nih.gov/28718724/) · DOI 10.1080/01635581.2017.1340491 |
| I007 | Sequential wheat, rice and oat staple periods | Human short sequential dietary study, 26 Mongolians | Rice lowered B. longum/B. adolescentis; wheat/oat favored some bifidobacteria; carbohydrate gene associations | Short, diet/population-specific; no general rice prescription to suppress B. longum | [PMID 28377764](https://pubmed.ncbi.nlm.nih.gov/28377764/) · DOI 10.3389/fmicb.2017.00484 |
| I008 | Commercial garlic powder, 0.1% and 1% w/v | In-vitro pure cultures and colonic model | Transient killing of B. longum DSMZ 20090/B. ovatus/C. nexile; L. casei DSMZ 20011 relatively resistant | Exact powder/strain and culture concentration; adaptation and recovery occurred; not all garlic or allicin formulations equivalent | [PMID 22480662](https://pubmed.ncbi.nlm.nih.gov/22480662/) · DOI 10.1016/j.phymed.2012.02.018 |
| I009 | Ferulic acid and gallic acid | In-vitro antimicrobial experiment | E. coli inhibition: ferulic MIC 100 micrograms/mL, gallic MIC 1,500 micrograms/mL; membrane effects | Keep gastrointestinal delivery/exposure uncertain; laboratory MIC is not a human oral dose | [PMID 23480526](https://pubmed.ncbi.nlm.nih.gov/23480526/) · DOI 10.1089/mdr.2012.0244 |
| I010 | Pistachio nuts | Rat streptozotocin-diabetes feeding experiment | Increased lactobacilli/bifidobacteria; reduced enterococci; changed community | Not a verified human E. coli eradication result; species-specific endpoint requires full tables | [PMID 32812934](https://pubmed.ncbi.nlm.nih.gov/32812934/) · DOI 10.1016/j.metop.2020.100040 |
| I011 | Edible vegetable blend oil versus lard | Mouse 42-day gavage experiment | Oil groups reduced E. coli, lactobacilli and bifidobacteria versus water control | Shows collateral changes and animal context, not selective pathogen suppression | [PMID 35034939](https://pubmed.ncbi.nlm.nih.gov/35034939/) · DOI 10.5650/jos.ess21247 |
| I012 | Garlic powder tablets, 400 mg containing 1,100 micrograms allicin/tablet, twice/day | Human double-blind trial; 43 randomized, 32 completed | Microbial changes were small/trends; Faecalibacterium/Bifidobacterium tended upward and Akkermansia downward | Both groups low-calorie diet; avoid presenting trend as proven genus increase | [PMID 36352899](https://pubmed.ncbi.nlm.nih.gov/36352899/) · DOI 10.3389/fnut.2022.1007506 |
| I013 | Low-glycinin or conventional soymilk versus bovine milk, 500 mL/day | Human randomized double-blind, 64 overweight/obese men | Bifidobacterium reduced in soy groups; overall taxonomic/diversity changes | Exact formulations differed; do not label soy universally Faecalibacterium-raising from abstract alone | [PMID 22895080](https://pubmed.ncbi.nlm.nih.gov/22895080/) · DOI 10.4161/gmic.21578 |
| I014 | Habitual chicken versus pork intake | Human observational comparison, 20 men/group | Different microbial profiles; chicken group had higher indole/skatole and SCFA | Not an intervention; higher indoles/BCFA not uniformly favorable | [PMID 34099832](https://pubmed.ncbi.nlm.nih.gov/34099832/) · DOI 10.1038/s41598-021-91429-3 |
| I015 | Cocoa husks replacing 7.5% of diet | Pig crossover feeding, 6 pigs | Increased F. prausnitzii; reduced Lactobacillus-Enterococcus and C. histolyticum groups | Cocoa husk feed is not interchangeable with chocolate/cocoa powder | [PMID 26877143](https://pubmed.ncbi.nlm.nih.gov/26877143/) · DOI 10.1021/acs.jafc.5b05732 |
| I016 | Genistein, daidzein, glycosides and equol | In-vitro study, 37 strains | Genistein/equol increased growth rate of F. prausnitzii and L. rhamnosus; equol inhibited some strains | Compound- and strain-specific dose effects; soy food cannot inherit every isolated-molecule effect | [PMID 28698467](https://pubmed.ncbi.nlm.nih.gov/28698467/) · DOI 10.3390/nu9070727 |
| I017 | Kiwi FFG standardized green kiwifruit powder | In-vitro SHIME using donor communities plus cell assays | Raised butyrate/F. prausnitzii/Roseburia/Bifidobacterium; AhR and cytokine effects | A human-derived gut model is not a human clinical trial | [PMID 37803696](https://pubmed.ncbi.nlm.nih.gov/37803696/) · DOI 10.1016/j.foodres.2023.113348 |
| I018 | Fresh ginger juice | Human short crossover study, 123 healthy adults | Community shifts; reduced Prevotella:Bacteroides and some Ruminococcus groups; sex differences | Do not turn any Fusobacterium increase into a health goal; taxa changes not uniform between sexes | [PMID 33708178](https://pubmed.ncbi.nlm.nih.gov/33708178/) · DOI 10.3389/fmicb.2020.576061 |
| I019 | Blueberry powder, 10% diet | Rat high-fat-diet experiment | Improved inflammation/insulin-signaling measures; increased Gammaproteobacteria | Useful polyphenol hypothesis; no human dose or desirable Fusobacterium target | [PMID 29490092](https://pubmed.ncbi.nlm.nih.gov/29490092/) · DOI 10.1093/jn/nxx027 |
| I020 | Habitual diet/habitat exposures | Human observational, 3,224 participants | Diet, geography, microbiome and serum metabolites associated | Cross-sectional/observational pathways do not establish individual food prescription | [PMID 38477427](https://pubmed.ncbi.nlm.nih.gov/38477427/) · DOI 10.1002/advs.202310068 |
| I021 | Walnut/pumpkin dietary pattern | Review built around an autism case narrative | Mechanistic discussion and reported long-term case improvement | Secondary evidence; not a controlled organism-modification trial; unassigned reference only | [PMID 37960217](https://pubmed.ncbi.nlm.nih.gov/37960217/) · DOI 10.3390/nu15214564 |
| I022 | Konjaku flour versus lotus-root-starch placebo | Human randomized double-blind, 69 obese adults | 5-week anthropometric/metabolic changes; increased Roseburia and reported R. inulinivorans | The study publication also classified C. perfringens favorably; curate taxon interpretation independently | [PMID 35300378](https://pubmed.ncbi.nlm.nih.gov/35300378/) · DOI 10.3389/fcimb.2022.771748 |
| I023 | Rhubarb extract standardized by rhein content | Human double-blind randomized trial, 30 days | Improved constipation; increased Roseburia/Agathobacter especially when initially low; fecal SCFA increased | Extract/laxative effect matters; do not promise engraftment or generalized microbiome cure | [PMID 36499011](https://pubmed.ncbi.nlm.nih.gov/36499011/) · DOI 10.3390/ijms232314685 |
| I024 | Commercial water-soluble annatto extract | In-vitro food microbiology | Inhibited S. thermophilus at 0.63% v/v and several Gram-positive organisms; no activity on tested Gram negatives/yeasts | Food-preservation assay; distinguish from tocotrienol extract and gut treatment | [PMID 12801012](https://pubmed.ncbi.nlm.nih.gov/12801012/) · DOI 10.4315/0362-028x-66.6.1074 |
| I025 | B. longum BB536, 4 billion CFU plus L. rhamnosus HN001, 1 billion CFU daily | Human preliminary trial, 20 healthy Italian adults | Combination associated with Blautia wexlerae increase and later Akkermansia increase/R. gnavus decrease | Combination effect cannot be assigned separately to either strain; no placebo comparison; strain detection persisted at 1 month; not permanent engraftment | [PMID 28487606](https://pubmed.ncbi.nlm.nih.gov/28487606/) · DOI 10.3748/wjg.v23.i15.2696 |
| I026 | Encapsulated/uncoated E. faecium, 5 million CFU/g feed | Broiler chicken experiment, 48 birds | Colonization and SCFA changes; barrier-gene downregulation also occurred | No human strain/product inference; exact strain genome and virulence/AMR identity required before translational card | [PMID 38761463](https://pubmed.ncbi.nlm.nih.gov/38761463/) · DOI 10.1016/j.psj.2024.103808 |
| I027 | B. animalis subsp. lactis BB-12 probiotic yogurt | Human randomized athlete study, 20 women, 8 weeks | Lower E. coli and increased selected genera; fatigue questionnaire endpoint | Do not generalize to any B. animalis or all yogurts; abstract's fatigue-score direction and conclusion warrant full-text check | [PMID 37374905](https://pubmed.ncbi.nlm.nih.gov/37374905/) · DOI 10.3390/microorganisms11061403 |
| I028 | Enterococcus durans EP1 | Mouse and cell/PBMC experiments | Increased fecal IgA/F. prausnitzii and altered inflammatory expression | Preclinical exact strain; not a general Enterococcus probiotic recommendation | [PMID 28239378](https://pubmed.ncbi.nlm.nih.gov/28239378/) · DOI 10.3389/fimmu.2017.00088 |
| I029 | Agave salmiana aguamiel concentrate and saponin-rich extract | Mouse high-fat-diet experiment | Increased Akkermansia with metabolic improvement | Not interchangeable with commercial agave syrup; formulation-mapping correction | [PMID 27678062](https://pubmed.ncbi.nlm.nih.gov/27678062/) · DOI 10.1038/srep34242 |
| I030 | Vitamin A intervention | Human uncontrolled pediatric ASD study, 64 followed / 20 microbiomes | Retinol and biomarker changes; lower Bifidobacterium; no significant autism-scale change | Not justification to give vitamin A to reduce beneficial bacteria; nutrient indication and study population differ | [PMID 28938872](https://pubmed.ncbi.nlm.nih.gov/28938872/) · DOI 10.1186/s12866-017-1096-1 |
| I031 | Xylitol, 1, 3 or 5 g/L culture concentration | In-vitro dynamic child-donor gut simulator | Dose-related Blautia/Anaerostipes/Roseburia and butyrate changes; cell-barrier effect | No human clinical dose; genus effect not demonstrated B. wexlerae-specific colonization | [PMID 38537869](https://pubmed.ncbi.nlm.nih.gov/38537869/) · DOI 10.1016/j.fct.2024.114605 |
| I032 | Iron or iron-containing micronutrient powders, 3 months | Human randomized placebo-controlled, 923 Bangladeshi infants | No significant overall microbiome effect in primary adjusted analysis; exploratory pathogen increases | Not supportive evidence for E. coli reduction; retain as null/limiting evidence, not positive antimicrobial option | [PMID 39367018](https://pubmed.ncbi.nlm.nih.gov/39367018/) · DOI 10.1038/s41467-024-53013-x |
| I033 | Oligonol litchi-derived polyphenol | Human randomized double-blind, 38 adults with NAFLD, 24 weeks | Steatosis improvement within active arm; increases in Akkermansia/Faecalibacterium and other taxa | Within-arm significance is not automatically significant treatment-control difference; exact product/condition | [PMID 35889878](https://pubmed.ncbi.nlm.nih.gov/35889878/) · DOI 10.3390/nu14142921 |
| I034 | Isoflavones with/without probiotic or FOS | Human randomized double-blind, 39 postmenopausal women | F. prausnitzii and other groups increased; FOS bifidogenic; equol-excretion-dependent changes | Combination arms and baseline equol phenotype matter; not all isoflavones/products interchangeable | [PMID 16317121](https://pubmed.ncbi.nlm.nih.gov/16317121/) · DOI 10.1093/jn/135.12.2786 |
| I035 | Standardized poplar-type propolis polyphenol extract | In-vitro simulated digestion/fermentation, 5 donor communities | SCFA increases and community changes after digestion; antioxidant activity declined | Digestion experiment adds exposure plausibility but remains preclinical; product specificity/allergy context | [PMID 35248845](https://pubmed.ncbi.nlm.nih.gov/35248845/) · DOI 10.1016/j.biopha.2022.112759 |
| I036 | Xylooligosaccharides, 2% diet | Gestational-diabetes mouse model, 60 mice | Increased Akkermansia, improved barrier markers/insulin resistance | Pregnant-mouse findings not a human pregnancy recommendation or individualized XOS prediction | [PMID 38426554](https://pubmed.ncbi.nlm.nih.gov/38426554/) · DOI 10.1039/d3fo04681h |
| I037 | High-molecular-weight barley beta-glucan, 3 g/day | Human randomized crossover | Bacteroides increased; changes depended on molecular weight, low-MW arms differed | Great formulation-specific fiber seed; abstract does not establish desired Bifidobacterium suppression | [PMID 26904005](https://pubmed.ncbi.nlm.nih.gov/26904005/) · DOI 10.3389/fmicb.2016.00129 |
| I038 | XOS, 0.04% or 0.40% diet | Cat feeding experiment, 24 cats | Blautia/Clostridium XI/Collinsella increased; Bifidobacterium decreased; dose effects | Different direction than common human bifidogenic findings; label species/model explicitly | [PMID 33138291](https://pubmed.ncbi.nlm.nih.gov/33138291/) · DOI 10.3390/molecules25215030 |
| I039 | Purified exopolysaccharides from intestinal Bifidobacterium strains | In-vitro fecal-slurry fermentation | SCFA and community changes; B. pseudocatenulatum EPS enriched F. prausnitzii in some donors; B. longum EPS differed | EPS is not live probiotic species and polymers/strains differ; selective E. coli reduction not established by abstract | [PMID 18539803](https://pubmed.ncbi.nlm.nih.gov/18539803/) · DOI 10.1128/aem.00325-08 |
| I040 | Partially hydrolyzed guar gum, 5 g up to 3 times/day | Human 20-volunteer sequential PAGODA trial | Increased acetate/butyrate/Faecalibacterium; decreased Roseburia/Blautia; most effects vanished after stopping | Display tradeoffs and transient response; no healthy-genus universality | [PMID 32354152](https://pubmed.ncbi.nlm.nih.gov/32354152/) · DOI 10.3390/nu12051257 |
| I041 | Arabinoxylan-rich or resistant-starch-rich diet | Pig feeding, 30 females, 3 weeks | Arabinoxylan strongest butyrate increase and F. prausnitzii/Roseburia shifts; distinct digestion sites | Specific substrate differences support separate modules; animal effect sizes not human predictions | [PMID 25327182](https://pubmed.ncbi.nlm.nih.gov/25327182/) · DOI 10.1017/s000711451400302x |
| I042 | Chicory forage rich in pectin | Pig and broiler feeding experiments | Fiber digestion and phylotypes differed by species and substrate | Not generic chicory-root inulin; forage and animal context important | [PMID 24341997](https://pubmed.ncbi.nlm.nih.gov/24341997/) · DOI 10.1186/2049-1891-4-50 |
| I043 | Inulin and XOS | Ex-vivo human stool and bacterial-isolate experiments | Identified direct inulin responders and indirect E. lenta/G. urolithinfaciens stimulation; XOS more bifidogenic | Excellent substrate/strain/crossfeeding evidence; no guarantee prebiotic response from abundance alone | [PMID 38097563](https://pubmed.ncbi.nlm.nih.gov/38097563/) · DOI 10.1038/s41467-023-43448-z |
| I044 | Genistein diet | Humanized-mouse breast-cancer experiment | Community changes and delayed tumor growth | Human-donor microbiota in mice remains animal evidence; no direct adult treatment efficacy | [PMID 29267377](https://pubmed.ncbi.nlm.nih.gov/29267377/) · DOI 10.1371/journal.pone.0189756 |
| I045 | Tannin-based supplement, 14 days | Human randomized placebo-controlled, 124 hospitalized COVID patients | No clinical improvement or significant overall microbiome shift; MIP-1alpha reduction | Null clinical primary context; not positive proof of targeted Bifidobacterium suppression | [PMID 36467850](https://pubmed.ncbi.nlm.nih.gov/36467850/) · DOI 10.1016/j.jff.2022.105356 |
| I046 | Tabebuia impetiginosa isolated anthraquinone-2-carboxylic acid and lapachol | In-vitro disk-diffusion assay | Strong C. paraputrificum inhibition for one compound; weak E. coli/C. perfringens inhibition; several bifidobacteria spared | Not proof Pink Trumpet Tree reduces B. longum; exact constituents and intestinal exposure unresolved | [PMID 15713033](https://pubmed.ncbi.nlm.nih.gov/15713033/) · DOI 10.1021/jf0486038 |
| I047 | Thyme/thymol and other spices | Review of food-preservation antimicrobial literature | Describes antimicrobial potential and delivery technologies | Secondary source; primary strain/MIC records required before effect claim; not a human gut trial | [PMID 37895922](https://pubmed.ncbi.nlm.nih.gov/37895922/) · DOI 10.3390/ph16101451 |
| I048 | Enzymolysis seaweed powder, 20 g/kg feed versus S. boulardii | Kitten feeding experiment, 30 kittens, 4 weeks | Faecalibacterium and barrier markers improved; SCFA not changed | Processed powder not generic seaweed; kitten-to-human translation uncertain | [PMID 37233678](https://pubmed.ncbi.nlm.nih.gov/37233678/) · DOI 10.3390/metabo13050637 |
| I049 | Laminarin, 300 ppm or fucoidan, 250 ppm-enriched seaweed extract | Pig feeding experiment, 75 pigs | Laminarin raised Faecalibacterium/Roseburia and lowered Campylobacter; both extracts increased butyrate | Product-specific animal data; lower alpha diversity despite some favorable outcomes | [PMID 33810463](https://pubmed.ncbi.nlm.nih.gov/33810463/) · DOI 10.3390/md19040183 |
| I050 | Plant-based diet, 12 weeks | Human 14-patient Crohn disease pilot | Faecalibacterium/Bacteroides increased; calprotectin decreased; stool plant DNA measured adherence | Small uncontrolled disease-specific study; not evidence for universal species restoration | [PMID 39545044](https://pubmed.ncbi.nlm.nih.gov/39545044/) · DOI 10.3389/fnut.2024.1502967 |
| I051 | Gluten-free diet then gluten/placebo challenge | Human study, 31 women with autoimmune thyroiditis | Bifidobacterium declined, Proteobacteria/other groups increased; investigators urged caution | Not an improvement card solely because bifidobacteria decrease | [PMID 38474814](https://pubmed.ncbi.nlm.nih.gov/38474814/) · DOI 10.3390/nu16050685 |
| I052 | Ketogenic diet | Rat ADHD model | Bifidobacterium/Blautia increased with behavioral and brain-neurochemical changes | No basis to tell adult to avoid KD just to suppress Bifidobacterium | [PMID 37585373](https://pubmed.ncbi.nlm.nih.gov/37585373/) · DOI 10.1371/journal.pone.0289133 |
| I053 | Gluten-free diet, 1 month | Human pilot, 10 healthy adults | Reduced B. longum/F. prausnitzii/lactobacilli; increased E. coli/Enterobacteriaceae | Useful counterevidence to simplistic diet advice; clinical celiac indication independent | [PMID 19445821](https://pubmed.ncbi.nlm.nih.gov/19445821/) · DOI 10.1017/s0007114509371767 |
| I054 | High-fiber dietary intervention | Human pediatric metagenomic reanalysis, 17 with PWS + 19 with simple obesity | F. prausnitzii strain SNP shifts and increased Bifidobacterium; nutrition-related gene changes | Community replacement/strain selection not proof directed mutagenesis; pediatric baseline-specific | [PMID 34135881](https://pubmed.ncbi.nlm.nih.gov/34135881/) · DOI 10.3389/fmicb.2021.683714 |
| I055 | Habitual dairy/plant intake | Human observational, 37 children aged 2–3 years | Food associations with diversity and taxa; apples/pears inversely associated with R. gnavus relatives | Not adult randomized dairy-avoidance evidence | [PMID 27694811](https://pubmed.ncbi.nlm.nih.gov/27694811/) · DOI 10.1038/srep32385 |
| I056 | Mediterranean versus low-fat high-carbohydrate diet, 2 years | Human dietary intervention, metabolic syndrome | Mediterranean arm partly restored F. prausnitzii, B. adolescentis/B. longum and other taxa | Long-term dietary pattern, not one food/gene-specific correction | [PMID 26376027](https://pubmed.ncbi.nlm.nih.gov/26376027/) · DOI 10.1016/j.jnutbio.2015.08.011 |
| I057 | Long-term moderate versus short high-intensity exercise | Mouse training and FMT experiment | Different taxa/SCFA and post-exhaustion muscle injury; transferability shown in mice | Not evidence all exercise raises Roseburia; do not extrapolate to post-exertional malaise patients | [PMID 39063080](https://pubmed.ncbi.nlm.nih.gov/39063080/) · DOI 10.3390/ijms25147837 |
| I058 | Chinese-herb complex during intermittent cold exposure | Female rat experiment | Complex improved barrier markers and increased Roseburia | Not human cold-exposure harm or proof avoiding cold corrects a stool taxon | [PMID 36532488](https://pubmed.ncbi.nlm.nih.gov/36532488/) · DOI 10.3389/fmicb.2022.1065780 |
| I059 | High-protein low-carbohydrate hypocaloric weight-loss program with/without live microorganisms | Human real-world program, 263 enrolled, 163 paired samples | At 10% weight loss, Akkermansia/P. distasonis increased; S. thermophilus/bifidobacteria/E. rectale decreased | Multi-component nonrandomized context; don't equate every decreased taxon with benefit | [PMID 35052696](https://pubmed.ncbi.nlm.nih.gov/35052696/) · DOI 10.3390/biomedicines10010016 |
| I060 | High versus low dairy, 6-week periods | Human randomized crossover, 46 overweight adults | High dairy increased S. thermophilus/E. ramosum; lower F. prausnitzii/B. wadsworthia; some constipation | Specific tradeoffs; predicted pathways did not change; no universal exclusion of fermented foods | [PMID 32708991](https://pubmed.ncbi.nlm.nih.gov/32708991/) · DOI 10.3390/nu12072129 |

### 17.1 Priority bindings and evidence corrections

- Fiber/substrate response cards: oats, PHGG, molecular-weight-specific barley beta-glucan, inulin/XOS crossfeeding, arabinoxylan/resistant starch, and kiwifruit powder.
- Glucosinolate mechanism/action card: broccoli sprouts trial, keeping plant myrosinase and microbial conversion distinct.
- Polyphenol conversion cards: walnuts with measured secondary bile acids; isoflavones/equol phenotype; standardized propolis and oligonol with model/product distinctions.
- Exact probiotic combination/strain cards: BB536 + HN001 combination, BB-12 yogurt, and preclinical EP1; never transfer one product result to a whole species.
- Antimicrobial candidates: garlic powder, ferulic/gallic acids and annatto extracts with specific in-vitro targets and concentration evidence, no invented oral dosing.
- Do not seed positive targeted-suppression cards from the iron, tannin, unspecified E. faecium, generic agave syrup, or cold-avoidance interpretations; the source studies do not establish those claims.
- No new vitamin B1/B3/B5/B6 deficiency treatment follows from microbial pathway percentiles; link gene-capacity explanations and separately measured host nutrient evidence.


## 18. Embedded synbiotic seed registry

Normalize these **16 exact strain–substrate records** and the whole-formulation record below into the existing intervention registry. The table preserves study exposures for reproducibility; these are **not automatically recommended personal doses**. Sources with limited extraction are marked explicitly. Complete their missing fields from the cited source before claiming exact formulation matching, while retaining supported evidence immediately. Human, animal, culture and model outcomes have separate evidence types.

| ID | Exact strain + substrate | Experiment / formulation | Finding and classification | Primary source / data |
|---|---|---|---|---|
| SYN001 | **Bifidobacterium adolescentis IVS-1 (iVS-1) + GOS** | Human double-blind six-arm RCT; 114 randomized, 94 analyzed, obese adults, 3 weeks. 10^9 CFU/day; 5 g/day GOS supplied as 6.9 g Vivinal powder (72.5% GOS, 22.8% lactose, 4.7% monosaccharides). Strain alone, GOS alone, combination and lactose placebo available. | Higher fecal abundance than BB-12; GOS did **not significantly enhance** strain abundance or clinical function. Some permeability markers improved, but **no functional synergism**. Critical negative-control benchmark despite successful earlier rat results. | https://doi.org/10.1186/s40168-018-0494-4 ; microbiome PRJNA434249 / SRP133159. |
| SYN002 | **B. animalis subsp. lactis BB-12 + same Vivinal GOS** | Same RCT and doses; independent BB-12 and BB-12+GOS arms. | No significant enhancement of fecal BB-12 with GOS; no functional synergism. Components showed some permeability effects, endotoxemia markers unchanged. Preserve as a separate strain-specific negative pair. | https://doi.org/10.1186/s40168-018-0494-4 . |
| SYN003 | **B. animalis subsp. lactis Bi-07 + XOS** | Human double-blind randomized factorial crossover; 44 recruited, 41 completed; 21 days per treatment. 10^9 CFU/day + 8 g/day xylo-oligosaccharides; probiotic, substrate, combination, maltodextrin control. | XOS bifidogenic, but fecal B. lactis not significantly higher in combination versus Bi-07 alone: **no specific ecological synergy**. Immune markers changed; these are not established clinical benefit. Carryover meant some immune outcomes used first period only. | https://doi.org/10.1017/S0007114513004261 ; https://pubmed.ncbi.nlm.nih.gov/24661576/ ; full text https://www.cambridge.org/core/journals/british-journal-of-nutrition/article/xylooligosaccharides-alone-or-in-synbiotic-combination-with-bifidobacterium-animalis-subsp-lactis-induce-bifidogenesis-and-modulate-markers-of-immune-function-in-healthy-adults-a-doubleblind-placebocontrolled-randomised-factorial-crossover-study/E0E19DDF4AAFA40CCB41E31E5CE670FA |
| SYN004 | **Lactobacillus acidophilus NCFM (ATCC 700396) + cellobiose** | Human randomized double-blind crossover, 18 healthy adults; 10^9 CFU + 5 g/day 97% cellobiose for 3 weeks, 3-week washout. Combination versus maltodextrin; **no component-only arms**. | Lactobacilli/bifidobacteria and **branched-chain** fatty acids increased; diversity and SCFAs did not. Do not label a demonstrated synergistic health effect or confuse BCFA with SCFA benefit. | https://doi.org/10.1111/1574-6941.12397 ; https://academic.oup.com/femsec/article/90/1/225/2680511 ; NCT01716910. |
| SYN005 | **Limosilactobacillus reuteri DSM 17938 + GOS** | 2014 human randomized crossover study compared GOS, rhamnose, combined substrates with probiotic. Same strain also tested with GOS formula in 82 Thai children aged 8–14 months: single-serving randomized isotope crossover; GOS 400 mg/100 mL formula, 235 mL serving, control and 2′FL-formula comparators. | 2014: substrate addition did not increase fecal counts or post-dosing persistence. **No added ecological benefit demonstrated**. Keep the later iron study as a separate endpoint/formulation, not evidence overturning persistence result. 2014 full numerical design should be re-extracted before production ingestion; primary abstract verified. | https://www.sciencedirect.com/science/article/pii/S1756464614001881 ; author thesis https://digitalcommons.unl.edu/foodscidiss/62/ ; 2024 https://pubmed.ncbi.nlm.nih.gov/39179207/ |
| SYN006 | **B. longum subsp. infantis EVC001 + purified LNnT** | SYNERGIE human single-blind randomized pilot: 62 Bangladeshi infants aged 2–6 months with severe acute malnutrition; 28 days 8×10^9 CFU/day alone (n20), same +1.6 g/day lacto-N-neotetraose (n21), lactose placebo (n21); follow-up day56. **No LNnT-only arm**. | EVC001 increased in both active arms, later declined. Probiotic improved WAZ and MUAC; combination improved WAZ versus placebo; **synbiotic superiority unproven**. Low breast-milk consumption and geography/diet mattered; colonization lower than healthy infants. Strong caution against extrapolating US breastfed-infant findings. | https://doi.org/10.1126/scitranslmed.abk1107 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC9516695/ ; ENA **PRJEB45396** (genomes, 16S, COPRO-seq, mouse RNA-seq). |
| SYN007 | **B. longum subsp. infantis LMG 11588 + six-HMO blend** | Ex vivo 48h colonic reactors; same six donors at ~3 months and ~12 months. Four core conditions blank / probiotic / substrate / combo; fifth adds B. lactis CNCM I-3446. 10^8 CFU/mL; 2.5 g total HMO/L. Mix: 3FL,2FL,DFL,3SL,6SL,LNT. Infant wt% 13.5,50.1,5.9,6.0,8.2,16.3; toddler wt% 37.4,29.3,3.4,14.4,5.3,10.2. | Combination increased HMO consumption and SCFA production, including initially low responders. **Ex vivo interaction**, no host health outcome. Substrate mixture/age must be encoded. | https://doi.org/10.1038/s42003-024-06628-1 ; **PRJEB71952**, ERX11887871–ERX11888086; processed sequencing Supplementary Data2. |
| SYN008 | **B. longum subsp. infantis M-63 + 2′FL** | 2026 ex vivo 24h pH-controlled YCFA cultures; four newborn donors (5–7d) and four children (3–5y); 1% w/v 2′FL (>91% purity), 10^8 cells/mL M-63. Control, 2′FL, 2′FL+M-63; **no M-63-only comparator**. | More bifidobacteria, acetate and aromatic lactic acids with combo than substrate; 2′FL alone not significant. Useful donor-context seed, but **cannot isolate formal synergy**; no clinical trial. Raw data on request only. | https://doi.org/10.3389/fnut.2026.1744839 |
| SYN009 | **B. longum subsp. infantis EFEL8008 + 2′FL** | 2025 in vitro digestion plus adult fecal fermentation; 19 volunteers supplied stool. Four conditions including both components and control; 3.51×10^6 CFU/mL +1% w/v 2′FL; 12/24h. | Combo increased bifidobacteria/SCFAs and reduced betaine-to-TMA conversion in this system. **No human intervention outcome**. qPCR targeted subspecies sialidase, not unique strain; technical triplicates do not establish 19 independent donor responses (aggregation requires checking). Raw data on request. | https://doi.org/10.20517/mrr.2025.35 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC12540058/ |
| SYN010 | **Anaerostipes caccae LAHUC + lactulose** | 2024 primary study in gnotobiotic allergic-infant-microbiota mice and antibiotic-treated SPF mice; component/control comparisons. Lactulose 5 g/L in drinking fluid; strain gavage 10^6 CFU (gnotobiotic), 10^7 CFU (SPF). | Luminal butyrate and allergy protection in prophylactic/therapeutic mouse models. **Engraftment alone did not increase butyrate**; a competing resident strain impeded engraftment. Cross-feeding matters. Preclinical LBP candidate, not human allergy treatment. LAHUC is not interchangeable with DSM14662 / HM220. | https://doi.org/10.1016/j.chom.2024.05.019 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC11239278/ ; LAHUC genome **PRJNA1102404**, 16S **PRJNA1091036**. |
| SYN011 | **Limosilactobacillus reuteri BNCC 186135 + GOS** | 2024 mouse DSS colitis study, n6/group: components alone, combination, disease/control arms; 2×10^8 CFU/day and 5% w/w dietary GOS for 5 weeks, DSS final week. | Combination improved colitis/barrier outcomes and enriched Bacteroides acidifaciens/pentadecanoic acid pathway. **Mouse host-response synergy candidate**, not human evidence. Strain identifier in methods; #1/#2 additional mouse isolates also studied, distinguish experiment-specific strains. | https://doi.org/10.1038/s41467-024-53144-1 ; **PRJNA1048359**, **PRJNA1154812**, downloadable source data. |
| SYN012 | **B. animalis subsp. lactis CNCM I-3446 + BMOS** | 2020 48h infant fecal cultures; substrate screening 10 donors, detailed component/combo comparison two donors; also strain pre-grown on dextrose versus BMOS. BMOS: demineralized-whey-derived mixture containing GOS and natural 3′/6′SL; not pure HMO. | BMOS promoted strain abundance; metabolic advantage strongest low-bifidobacteria donor. Acetate often similar to BMOS alone; gas increased. **Ex vivo formulation/priming interaction**, not clinical efficacy. Exact composition/dose needs extraction for production use. | https://doi.org/10.3390/nu12082268 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC7468906/ |
| SYN013 | **L. acidophilus NCFM + enzymatically hydrolyzed oat β-glucan (OBGH)** | 2012 four-stage human-colon culture model after carbohydrate growth screen. | Growth and SCFA effects; exploratory **in vitro** candidate. Hydrolysate enzyme/linkage/DP matter; cannot substitute intact oat β-glucan. | https://doi.org/10.1371/journal.pone.0047212 |
| SYN014 | **B. animalis subsp. lactis Bl-04 + xylobiose** | Same 2012 model/screen. | Increased Bl-04 in culture and altered metabolites; not established host benefit. | https://doi.org/10.1371/journal.pone.0047212 . |
| SYN015 | **B. animalis subsp. lactis Bl-04 + raffinose** | Same 2012 model/screen. | Candidate for substrate–strain compatibility screening. Original model did not establish human efficacy or a clinical factorial interaction. | https://doi.org/10.1371/journal.pone.0047212 . |
| SYN016 | **B. longum subsp. infantis PLMB0001 + donor-milk-derived HMO concentrate** | 2023 **unblinded** healthy-adult antibiotic-perturbation study, 56 participants, three cohorts: antibiotics; antibiotics+probiotic; antibiotics+probiotic+HMO. Commercial probiotic isolate explicitly designated PLMB0001; ≥8×10^9 CFU/day days1–14 and 18 g/day HMO days1–28. No HMO-only arm. | High, variable substrate-dependent engraftment and metabolite/community changes; **no clinical efficacy endpoint established**. Crucial identity correction: PBI001 was a separate isolate used in mouse/co-culture experiments, not the adult administered commercial probiotic. | https://doi.org/10.1016/j.chom.2023.08.004 ; https://pubmed.ncbi.nlm.nih.gov/37657443/ ; accepted fulltext https://escholarship.org/content/qt3p87s1jc/qt3p87s1jc.pdf ; **PRJNA993161**, NCT05141903. |


### 18.1 Whole-formulation record and identity controls

**SYN017 — DS-01 bundle [S14]:** studied preparation contains 53.6 billion AFU (not CFU), 400 mg pomegranate extract (>40% polyphenols) and 24 strain labels below. Positive gastrointestinal outcomes apply to the tested bundle; no component arms, personalized allocation or microbiome measurement establishes which component caused them. Ingest the original section-2.3 identities and exposure protocol; keep the bundle atomic.

| Publication taxon | Strain labels in the tested formulation |
|---|---|
| B. longum | SD-BB536-JP; HRVD90b-US; SD-CECT7347-SP |
| B. breve | SD-BR3-IT; HRVD521-US |
| L. plantarum | SD-LP1-IT; SD-LPLDL-UK |
| L. rhamnosus | SD-LR6-IT; HRVD113-US; SD-GG-BE |
| B. infantis | SD-M63-JP |
| B. lactis | SD-BS5-IT; HRVD524-US; SD150-BE; SD-CECT8145-SP; SD-MB2409-IT |
| L. crispatus | SD-LCR01-IT |
| L. casei | HRVD300-US; SD-CECT9104-SP |
| L. fermentum | SD-LF8-IT |
| L. reuteri | RD830-FR; SD-LRE2-IT |
| L. salivarius | SD-LS1-IT |
| B. adolescentis | SD-BA5-IT |

Preserve original names alongside resolved taxonomy. Supplier-style labels require independently verified aliases before being merged with culture/genome identities or with other products. Do not infer component quantities where only a total is published.

**Genome anchors for the seed registry:**

- iVS-1: **CP123050**; BioProject **PRJNA954850**, BioSample **SAMN34156685**, raw **SRR24147545**. https://doi.org/10.1128/mra.00541-23 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC10720497/
- BB-12: **CP001853**. https://doi.org/10.1128/JB.00109-10 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC2863482/
- Bl-04: **CP001515**, distinguished from DSM10140 **CP001606**. https://doi.org/10.1128/JB.00155-09 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC2698493/
- NCFM old reference **NC_006814** supported by primary molecular papers, but 2024 re-curation found differences from historical sequence; prefer lot/assembly validation. https://pmc.ncbi.nlm.nih.gov/articles/PMC11697832/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC7858073/


These accession roots must be resolved to exact current versions and downloaded sequence hashes during ingestion. A historic reference genome is not proof that every commercial lot contains that exact genome.

**Required identity/evidence corrections:** do not fabricate AH1206+GOS or 35624+fibre trials from brand associations or probiotic-alone studies. A fermented carrier is not automatically a tested selective substrate. PLMB0001, PBI001, EVC001 and ATCC15697 are not interchangeable strain IDs; high ANI alone does not make them identical. LAHUC is not DSM14662/L1-92. A study of breastfed infants is not a purified-HMO factorial trial. The 2012 growth screen's digestible sugars do not automatically qualify as prebiotics.

## 19. Targeted intervention discovery across all legitimate evidence

### 19.1 Binding inclusion and report behavior

The recommendation engine must search and use **controlled laboratory experiments, biochemical assays, isolated strains, co-cultures, human-derived ex-vivo communities, organoids, animal experiments, human observations and human trials**. A human trial, regulatory approval, commercial product, durable engraftment or direct measurement of intestinal concentration is **not** required for a supported candidate to appear. Evidence setting determines the claim that can be made; it does not decide whether the finding is worth showing. Include a relevant experimental compound even when the tested delivery system is a research formulation.

One relevant, traceable favorable experiment is enough to create a candidate card. Preserve methodological quality, replication, null findings and opposing results alongside it. A biochemical mechanism, genome prediction or simulation alone creates a **mechanistic candidate**, visibly distinguished from observed inhibition, growth, metabolism or host improvement. A retailer claim without a traceable experiment is a discovery lead; it cannot establish the effect. A user's observed response can be stored as a supplied personal observation without being converted into a controlled study or discarded.

The core output is **targeted options with evidence**, not a claim that every detected organism should be eliminated. Options from nonhuman evidence belong in the normal finding-specific option list with their evidence labels. They must not be hidden by a default human-only filter or placed in a collapsed section solely because no human trial exists. A broader evidence library remains searchable. Confirmed contraindications, preparation incompatibility, unsupported target mapping and documented interactions have explicit, independently sourced effects on prioritization; absence of human research is not a contraindication.

Do not emit “untreatable,” “nothing can help,” “no effective supplements,” or “no options” merely because a human trial was not found. The precise alternatives are `search_not_run`, `search_incomplete`, `no_relevant_record_found`, `supported_options_found`, `mechanistic_candidates_only`, or `evidence_conflicted`. A complete search of a specified snapshot is not a search of every experiment ever performed.

Use all existing finding types as retrieval triggers: present/high/low/unusual organisms, resolved strains/pathotypes, pathogen determinants, function capacities, biofilm components, mycobiome findings, existing disease-pattern outputs, symptoms and supplied external measurements. Preserve each finding's actual measurement and interpretation. Do not change any existing score to manufacture a recommendation.

### 19.2 Required finding-to-option routes

| Trigger | Mandatory search/retrieval routes | Required distinction on the card |
|---|---|---|
| Concerning organism expansion or pathogen signal | Target-specific compounds and botanical preparations; antagonistic probiotic strains; colonization-resistance communities; substrate/environment adjustments; applicable clinician-review treatments | Species versus tested isolate/pathotype; growth inhibition versus killing, prevention, carriage reduction or symptom improvement; relative abundance versus absolute burden |
| Biofilm-associated findings | Target- and matrix-specific attachment prevention, quorum/adhesion effects, matrix disruption, viable-cell reduction, dispersal, and combination studies | Existing DNA score is a signal of genetic potential; it does not measure a physical mature biofilm. Reduced crystal-violet biomass is not necessarily killing; detachment is not necessarily clearance |
| Low supported butyrate or other desirable capacity | Exact live producers, cross-feeding partners, substrates for resident producers, measured consortia, relevant postbiotic or metabolite-delivery alternatives | Genomic capacity, culture output, fecal concentration and host benefit are different endpoints; stool concentration is influenced by production, absorption and transit |
| Low or undetected potentially useful organism | Exact-strain preparations if available, experimentally supported substrates, co-culture/cross-feeding options and functionally substitutable routes | Non-detection is not proof of absence or deficiency; a commercial strain need not permanently engraft to have an effect |
| Barrier-associated microbial findings or relevant supplied symptoms/tests | Zinc carnosine/polaprezinc, preparation-specific glutamine studies, relevant strains/substrates, lactoferrin and experimentally supported peptides | Separate actual host permeability experiments from stool-DNA proxies. Preserve oral, injected, topical and engineered delivery routes |
| High fungal/opportunistic fungal signal | Exact fungal-target antifungal and antibiofilm assays, antagonistic organisms, preparation-specific host/community studies | Gut colonization is not invasive infection; a vaginal/oral assay cannot silently become an intestinal treatment result |
| Disease-pattern or symptom goal | Condition-specific host effects plus separately supported microbial mechanisms | Similarity is not diagnosis; targeting an associated organism does not establish disease modification |
| Existing antimicrobial exposure or resistance finding | Susceptibility-qualified experimental options, complementary probiotic evidence, interactions and ecological recovery | Resistance genotype cannot establish sensitivity to a botanical; an unlinked AMR gene is not a resolved resistant strain |

These routes are additive. A card may link to several findings, but source effects must not be counted repeatedly as independent evidence. Show useful alternatives when one ingredient has an unfavorable target-specific tradeoff. Do not automatically assemble every individually promising option into a simultaneous regimen.

### 19.3 Federated research assets and ingestion

No single database establishes every probiotic's identity, delivery, antimicrobial activity and condition-specific benefit. Implement a federated registry using the resources and primary seeds embedded here. Human-indication guides are one layer; natural-product assays, strain traits, growth matrices and experimental papers are equally necessary input layers.

| Resource | Official source / primary record | Required adapter behavior |
|---|---|---|
| NPASS 3.0 / 2026 | [Primary database paper](https://pubmed.ncbi.nlm.nih.gov/41243954/); [database](https://bidd.group/NPASS/index.php); [downloads](https://bidd.group/NPASS/downloadnpass.html) | Import compound identity, source organisms, experimental target/activity, preparation, references and separate ADME/toxicity records. The publication describes 204,023 natural products and 1,048,756 activity records across biological applications; these are not counts of gut treatments. Pin the actual release and verify its terms. |
| PubChem BioAssay and PUG-REST | [BioAssay documentation](https://pubchem.ncbi.nlm.nih.gov/docs/bioassays); [PUG-REST](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest); [tutorial](https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest-tutorial) | Preserve AID, SID and CID separately, protocol, concentration/unit, active/inactive/inconclusive result, endpoint and source publication. “Active” alone is insufficient for a gut-target claim. Import counter-screens and inactive results. |
| ChEMBL | [API documentation](https://www.ebi.ac.uk/chembl/api/data/docs); [web-services usage](https://github.com/chembl/GLaDOS-docs/blob/master/web-services/chembl-data-web-services.md); [assay/target semantics](https://github.com/chembl/GLaDOS-docs/blob/master/frequently-asked-questions/chembl-data-questions.md) | Import source, document, assay, target organism, molecule/form, activity relation/value/unit and quality flags. Its confidence score concerns **target assignment**, not therapeutic efficacy. Whole-organism assays can have score 1; do not discard them using a protein-target threshold such as ≥8. |
| B-AMP | [Official resource and terms](https://b-amp.karishmakaushiklab.com/); [v2 primary paper](https://doi.org/10.3389/fcimb.2022.1020391) | Discovery of antimicrobial peptides, biofilm targets and structural models. Distinguish empirical effects from docking. A modeled interaction is not a demonstrated antibiofilm effect or an available oral product. Preserve separate data/code licenses. |
| Primary literature | [PubMed](https://pubmed.ncbi.nlm.nih.gov/); [Europe PMC REST documentation](https://europepmc.org/RestfulWebService); publisher supplements | Search actual target and ingredient aliases across all experimental settings, with no human-only/RCT inclusion filter. Store source document, correction status, extraction location and exact supported claim. Do not redistribute restricted full texts. |
| Probiotic, strain, utilization and product resources | Resource register and exact-strain seed records below; S01–S15; SYN001–SYN017 | Join by verified strain/preparation identities. Search identity, delivery, antagonism, substrate utilization, cross-feeding and host effects independently rather than expecting one trial to provide every axis. |

Deduplicate imported experiments across NPASS, ChEMBL, PubChem, reviews and the original paper. Use DOI/PMID plus experimental arm, preparation, target isolate, endpoint and observation identifier; an assay replicated across databases is one experiment. Keep database record IDs as provenance aliases. Verify chemicals by structure/InChIKey and stereochemistry where relevant; do not merge salts, pure compounds, standardized extracts and whole herbs into one efficacy entity.

Implement resumable adapters with raw-source hashes, pagination, schema/version checks, polite request limits, error states and local replay fixtures. A failed online lookup is not a negative literature result. Cached reports must work without internet. Source terms and restricted APIs cannot become a hidden mandatory dependency; manually curated source-backed records remain usable. A database absence does not mean no evidence exists.

### 19.4 One observation, several independent evidence axes

Extend the section-8 evidence edge rather than building a second disconnected recommendation system. The following fields are required in the normalized observation and candidate join:

- **Identity:** intervention entity, exact chemical/strain/preparation, mixture constituent identity, target organism/isolate/pathotype, publication taxonomy, current taxonomy, source body site and tested community. Identity resolution includes a reason, not just a numerical similarity cutoff.
- **Experiment:** setting, design, controls, independent biological units, technical replicates, assay medium/oxygen/pH, duration, measured endpoint, direction, uncertainty, source locator and extraction status. Store unreported values as null. `biochemical`, `in_vitro`, `ex_vivo`, `organoid`, `animal`, `human_observation`, `human_trial` and `in_silico` are distinct allowed settings.
- **Effect endpoint:** growth support/inhibition, CFU or viability, bactericidal/fungicidal effect, mature biomass, attachment prevention, dispersal, metabolite generation, absolute/relative abundance, carriage persistence, permeability, symptoms or disease endpoint. Preserve reported MIC/MBC/MBIC/MBEC definitions and units; labels are not interchangeable.
- **Delivery:** tested route, intended site, acid/bile survival, dissolution/release evidence, simulated digestion, live human recovery, concentration at site when measured, formulation rationale, absorption/metabolism and unresolved exposure. Local gut action and systemic bioavailability are separate.
- **Persistence:** detection during dosing, viable recovery during dosing, persistence after washout, duration, baseline presence and exact-strain resolution. Neither persistence nor durable engraftment is mandatory for candidate inclusion.
- **Tradeoffs:** target-matched stimulation of concerning organisms, measured effects on beneficial organisms, tolerance, allergens, toxicity and interaction evidence, contradictory results and assessment coverage. Unknown remains unknown, not harmful and not safe by default.
- **Translation:** exact measured statement, proposed relevance to the user's finding, target/route/population mismatch, availability as tested, commercial preparation bridge and what remains uncertain. Human benefit cannot be inferred from a label or genomic trait alone.

The recommendation must expose these axes separately. An organism can survive transit without a demonstrated desired effect; a metabolite or inactivated preparation can act without being alive; an ingredient may inhibit a target in culture with uncertain intestinal exposure. All three are valid distinct records.

A delivery evidence state is one of `human_site_measurement`, `human_live_recovery`, `animal_site_measurement`, `direct_site_administration`, `oral_outcome_observed`, `strain_dna_recovery_only`, `simulated_digestion`, `acid_bile_assay`, `formulation_rationale`, `not_applicable_nonviable`, `route_mismatch`, `unknown`, or `contradicted_for_tested_preparation`. More than one observation may support or conflict with a candidate. Do not reduce them to one unqualified “reaches the gut” flag. A species-level stool sequence during dosing is insufficient to establish live recovery of the commercial strain. Normalize seed wording at import: `viable_fecal_recovery` maps to `human_live_recovery` only for the human culture evidence recorded; `simulated_GI_survival` maps to `simulated_digestion`; `acid_bile_screen_only` maps to `acid_bile_assay`; `not_tested` maps to `unknown`. Preserve raw wording. An `engineered_enteric_release` label requires inspection of the actual assay before mapping to simulated digestion, site measurement or formulation rationale. `oral_outcome_observed` indicates route and outcome, not directly measured intestinal exposure. `not_applicable_nonviable` concerns bacterial viability only, not delivery of the active material.

Enteric or colon-targeted capsules are legitimate delivery hypotheses and research preparations. Link exact dissolution/release experiments where available. “Could be enteric-coated” creates `formulation_rationale`, not verified survival or a product match. Conversely, absent delivery studies do not erase experimentally observed activity: show **activity supported; intestinal exposure unresolved**. Low systemic absorption is not proof that unchanged active compound reaches the colon at an inhibitory concentration.

### 19.5 Relevance-first ranking across evidence settings

Apply `targeted-option-ranking/2` to newly added intervention options and section-11 synbiotic candidates. Remove blanket ordinal rules in which every human study outranks every animal or culture study. Retain setting prominently and match the study endpoint to the actual claim. A human constipation study cannot support eradication of Klebsiella; a direct Klebsiella inhibition experiment can support an experimental Klebsiella-targeted option.

For each candidate-goal edge, store the following **ordering codes**, which are not probabilities or biological effect sizes:

| Field | Values |
|---|---|
| `target_match` | 3 exact relevant resolved target/isolate/context; 2 relevant species/function/condition with explicit unresolved subtype; 1 related target or explicit mechanistic bridge; 0 unrelated/no supported bridge |
| `endpoint_match` | 2 actual desired endpoint observed; 1 measured intermediate relevant to the goal; 0 mechanism/model only |
| `preparation_match` | 3 exact studied formulation; 2 exact constituent(s) in a disclosed new formulation; 1 ingredient/species/preparation proxy; 0 unresolved |
| `empirical_support` | 2 independently reproduced relevant observation; 1 one qualifying experimental observation; 0 prediction/rationale only; null extraction incomplete |
| `delivery_match` | 2 relevant site/route evidence; 1 preparation-supported plausibility; 0 unknown or mismatched route, with these two states separately displayed |

`target_match=3` requires that the **sample** resolution supports the asserted match; a strain-specific paper and a species-only sample normally yield 2, not 3. A category-level target such as gut barrier must identify the intermediate or host endpoint and must not claim measured barrier injury from sequence data. “Exact endpoint” is judged at the claimed evidence setting: a culture experiment is exact for growth inhibition, not for human infection clearance. Target and endpoint matching are curator-controlled, versioned joins, not free-text semantic guesses.

For a multi-goal candidate, compare its vectors goal by goal; unsupported goals receive target/endpoint/support 0 and remain visibly uncovered. Use section-11 goal-coverage C and pair-coverage P for summaries only; neither is a dominance dimension or tie-break. P is null for single-component candidates. Do not collapse the best edge across different targets to create a falsely strong composite edge. Do not promote a candidate by adding weak duplicate studies. Independent replication requires separate biological experiments; shared-cohort papers and repeated technical wells do not qualify.

Use one ordered goal list G for the comparison scope. For each candidate and goal, retain all edges and choose one **supporting edge tuple** lexicographically by descending target, endpoint, preparation, empirical and delivery codes, with source/observation ID as final tie-break; a null sorts after every defined value for this choice. Every tuple dimension must describe that same evidence path; do not combine an unrelated clinical endpoint and an unrelated delivery study into a fictitious experiment. Compatible separately sourced delivery observations may attach through a verified preparation/route join and must retain their own provenance. Unsupported goals have an explicit all-zero tuple; unextracted observations have nulls, not fabricated zeros.

Compute nondominated fronts over the concatenated goal-specific five-code tuples. Two candidates compare only within the same goal list, eligibility class and evidence/model scope. Dominance requires all components defined, no smaller component, at least one larger component, and identical tradeoff-assessment coverage maps. In addition, the dominating candidate's sets of supported adverse endpoints and contrary observations must be subsets of the other candidate's sets. Different assessment coverage or null components makes that comparison incomparable. Same vectors are nondominating ties. Remove each nondominated front and repeat.

For deterministic ordering within a front, use supplied preference (preferred, neutral/unknown, avoid), reported intolerance (none/unknown before reported, with unknown labelled), then five descending user-weighted mean codes in this order: target, endpoint, preparation, empirical support, delivery. For dimension k, use `sum(w[j] * code[j,k]) / sum(w[j])` over the complete G, with default equal weights, and reject negative/nonfinite weights. A null in any positively weighted code gives a null mean; sort null means last using an internal ordering sentinel that is never saved or displayed as a measured value. Finish with ascending component count and lexicographic stable candidate ID. If all goal weights are zero, leave candidates unranked. Single ingredients and combinations use the same five dimensions; P is not needed for their comparison. Publish every join, mean and tie-break as an engineering relevance explanation, not an efficacy estimate.

Nonhuman setting and unknown delivery do not make a candidate ineligible or automatically hide it. Use explicit reasoned classes: `candidate_option`, `mechanistic_candidate`, `context_caution`, and `contraindicated_for_person`. Unknown commercial availability is an availability state, not lack of evidence. A research-only peptide or organism retains its finding-specific evidence card and can be marked **investigational; route/product unresolved** rather than receiving a personal self-use instruction. Confirmed contraindications are excluded from the start-here selection but remain visible with evidence. A user selecting a “human studies” filter may narrow the list; that must not be the default.

Safety/contraindication screening here is specific to the known ingredient, formulation and supplied patient facts. It is not a requirement that every candidate have pharmaceutical approval. The engine must not auto-generate human doses from MICs, animal experiments or injections, combine incompatible routes, or manufacture a claim that an experiment established cure.

### 19.6 Competitive feeding, selective effects and combination logic

Evaluate both **support desired organisms/functions** and **avoid supported adverse ecological effects**. Use exact substrate chemistry and experimentally matched strains where possible. HMOs, fibre, amino acids and lactoferrin are not uniformly beneficial or uniformly harmful to all communities. Do not treat their generic marketing category as selectivity evidence.

A concerning taxon's genome encoding utilization creates a hypothesis; cultured utilization creates a measured feeding edge; expansion in a host community adds a context-specific ecological observation. These remain separate. Uptake, host absorption, competing organisms, substrate accessibility, dose and existing diet can change the outcome. Do not conclude that oral glutamine feeds an undesirable organism simply because bacteria use nitrogen or carry a glutamine-related enzyme. Equally, do not suppress actual adverse feeding data when it exists.

For every proposed strain/substrate/compound combination check: exact supported target edges; exact-strain and ingredient overlap; antagonism between its components; overlapping adverse effects or allergens; available utilization/cross-feeding evidence; and whether the combination was actually tested. Unknown interactions remain visible. Do not infer additivity or synergy from two individually favorable papers. A direct antimicrobial and a live probiotic can conflict; show measured compatibility if available and otherwise unresolved interaction, without inventing administration timing.

Produce useful alternate paths when a feeding tradeoff exists: another supported substrate, a different exact strain, a host-directed option, a functional substitute or a researched postbiotic. The ranking must retain multiple paths rather than force one universal “best supplement.”

### 19.7 Retrieval, coverage and implementation outputs

For every rendered finding create an `option_lookup` record containing finding ID, canonical targets and uncertainty, registry/search snapshot, search strategy IDs, eligible candidate IDs, mechanistic candidate IDs, adverse/contrary edge IDs, unresolved identities, lookup state and last reviewed date. Coverage includes normal findings when the person explicitly selects a relevant goal; normality alone does not require intervention.

Required discovery query families are target synonyms × intervention classes × endpoint terms. Include taxonomic renamings, strain aliases, plant Latin names, chemical names/salts and preparation synonyms. Endpoint families include inhibit/growth/MIC, kill/viability, adhesion/biofilm/matrix, colonization/competition, metabolite/cross-feeding, permeability/tight-junction and host outcome terms. Run compound-first reverse searches as well as target-first searches. Search delivery and contrary/null outcomes separately. Do not add clinical or human filters to the universal query. Record searched source/version and coverage, rather than declaring exhaustive global knowledge.

Automated extraction produces candidate records with provenance and explicit unresolved fields. It must not invent missing strain IDs, doses or effect directions. Claims become available at their supported extraction level: a complete abstract can support its explicit target/model/result; unobserved numerical details require full-text or supplement verification. A source with credible contrary evidence remains indexed. Retractions/corrections change the affected claim status and all dependent views.

Expose these operations through the existing CLI conventions: refresh/pin intervention assets; import/validate evidence; explain one finding's candidate matches; list coverage gaps; render from cached records. Reuse existing commands where possible rather than imposing new root directories. Persist canonical registry, candidate and coverage objects as application outputs with schema versions and hashes. The implementation may create its normal code/data/test files; this specification is the sole handoff document.

### 19.8 Report contract and acceptance scope

Within **What you can do**, show one to three clearly prioritized options and a complete searchable list with the number of additional options. Each finding detail must link directly to its options. Do not stop at the top three or silently drop preclinical candidates. Integrate deeper substrate/strain reasoning into **Microbial Functions & Metabolism** and the existing relevant pathogen, organism, biofilm, mycobiome or condition detail page.

Each compact card answers: **what it targets; why it matches this finding; what was demonstrated; in what experimental setting; exact ingredient/strain/preparation; gut-delivery evidence; important contrary results/tradeoffs; and the primary source**. Use “inhibited growth in culture,” “reduced colonization in mice,” “improved the measured human endpoint,” or another exact statement. Never convert all of them to “proven to fix your gut.”

The resource and seed tables following this contract are mandatory initial coverage, not the ceiling of the registry. Ingest their supported facts, expose unresolved fields precisely and enable further targeted imports using the same schemas. Completion requires positive, null, conflicting, preparation-mismatch, delivery-unknown and exact-strain fixtures across several settings. Broad evidence inclusion and truthful endpoint labeling must pass together.


## 20. Exact-strain probiotic resources and initial target records

These records extend the existing I and SYN registries. Human studies and experimental models are independent evidence layers. Resource trait flags are discovery inputs; preserve exact materials, endpoints and unresolved fields. Records TP01–TP18 are mandatory initial coverage.

### 20.1 Probiotic resource adapters

| Resource | Useful content / implementation | Access, export, license and limits | Official/primary URL |
|---|---|---|---|
| Probio-Ichnos | Best immediate open structured seed for experimentally tested acid/bile resistance, adhesion, antimicrobial, immunomodulatory, antioxidant and antiproliferative properties at strain level. Records have PMIDs; distinguish present, absent, missing, uncertain. | Public GitHub includes `probioIchnos.json` and MIT LICENSE. Article reports 12,993 entries/11,202 strains/470 species, but repository README gives 12,887/11,074/443: preserve source snapshot/count rather than claiming these identical. Binary author-conclusion properties require primary-paper enrichment for exact pathogen, concentration, exposure, assay, null findings. No universal human-efficacy field. | https://github.com/Mtsif/probio-ichnos ; https://github.com/Mtsif/probio-ichnos/blob/main/probioIchnos.json ; https://github.com/Mtsif/probio-ichnos/blob/main/LICENSE ; paper https://doi.org/10.3390/microorganisms12101955 ; PMID 39458265 |
| IPDB (Integrated Probiotic Database), 2024 extension | Genomes and manually curated strain functional annotations; useful for adhesins, glycosidases, bacteriocin-production potential, metabolic capacity, product-label-to-strain discovery. 2024 adds 18 lactobacilli/Streptococcus/Heyndrickxia strains to 34 bifidobacterial strains. | Authors explicitly provide downloadable ZIP, freely accessible. No documented public API found. Specific dataset redistribution license not verified; paper license alone must not be assumed to cover every genome/source. Genetic capacity is `in_silico_predicted`, never measured antagonism or human benefit. Strain aliases require reconciliation (paper itself flags extremely similar genomes deposited under different labels). ZIP endpoint documented by paper, live binary retrieval not tested here. | http://probiogenomics.unipr.it/cmu/ ; http://probiogenomics.unipr.it/files/IPDB_latest.zip ; https://www.oaepublish.com/articles/mrr.2024.11 ; DOI 10.20517/mrr.2024.11 |
| BacDive (DSMZ) | Strain identifiers, culture collection crossrefs, taxonomy, physiology, metabolism, source, cultivation, sequence links, biosafety; useful ID backbone and supporting survival/metabolism facts. | Official REST API https://api.bacdive.dsmz.de/ and SPARQL https://sparql.dsmz.de/api/bacdive ; downloadable selections. About page reports CC BY 4.0 and requests commercial users contact DSMZ; also mentions DSI benefit-sharing. API v2 official client says no registration now required. Do not confuse MIT client license with data license. Not an intervention efficacy database. | https://www.bacdive.dsmz.de/about ; https://github.com/JKoblitz/bacdive-api ; https://sparql.dsmz.de/bacdive/ |
| ProBioQuest | Literature discovery/co-occurrence engine joining PubMed, ClinicalTrials.gov, PatentsView; supports exact strain keyword combinations and all evidence types. | Paper documents CSV export of search results and source link-outs; not a validated intervention graph. It retrieves source abstracts/summaries and infers relevance, not causal benefit. Public API or full database dump not verified. Article CC BY-NC 4.0; do not infer database/commercial redistribution rights from free search access. Current service operation not independently confirmed; use direct upstream APIs as resilient fallback. | http://kwanlab.bio.cuhk.edu.hk/PBQ/ ; https://academic.oup.com/database/article/doi/10.1093/database/baac059/6645125 ; PMID 35849028 |
| AEProbio USA Clinical Guide | Practical condition → studied exact strain/formulation → country-specific commercial product layer; human clinical evidence only, one layer among several. | Public website/app and PDFs; no public structured API or open bulk-reuse license found. Independently extract primary studies and verified product labels; do not treat unlisted preclinical candidates as unsupported. Current guide referenced by ISAPP March 2026. | https://usprobioticguide.com/ ; https://isappscience.org/resource/clinical-guide-to-probiotic-products-available-in-usa/ |
| AEProbio Canadian Clinical Guide | As above, Canadian product catalog and evidence grading. | https://aeprobio.com/resources/Canadian_Guide_En.pdf resolves to ~40.5 MB according to web retrieval, too large for extraction in this pass. No open API/license verified. Website may error intermittently; not grounds to invent a structured feed. | https://www.probioticchart.ca/ ; https://aeprobio.com/resources/Canadian_Guide_En.pdf |
| UK Probiotic Guide | Historical UK product mappings may aid leads, but **not currently maintained**. | Live homepage says guide no longer available; launched 2025-03-01, now outdated/on hold due to lack of ongoing funding. It points users to USA/Canada guides. Must date any historical products and verify current manufacturer identity. | https://probioticguide.uk/ |
| WGO Global Guideline 2023 | Tables of condition-specific strains/formulations and clinical evidence; useful reference extraction and conflicts with other guideline recommendations. | Public HTML/PDF; no public API or open data redistribution license verified. Latest probiotic guideline on WGO index is 2023. It is not a comprehensive preclinical database. | https://www.worldgastroenterology.org/guidelines/probiotics-and-prebiotics/probiotics-and-prebiotics-english ; https://www.worldgastroenterology.org/UserFiles/file/guidelines/probiotics-and-prebiotics-english-2023.pdf |
| NIH ODS professional fact sheet | High-level interpretation, strain specificity, product quality/viability issues, trial references, safety/context. | Free reference page; no strain-target structured database. Trace claims to trials. | https://ods.od.nih.gov/factsheets/Probiotics-HealthProfessional/ |
| ClinicalTrials.gov API v2 | Registered human interventions, exact strain/formulation mentions, comparators, primary/secondary endpoints, enrollment, status, results if submitted, linked publications; null/unpublished evidence discovery. | Official JSON REST, CSV downloads and OpenAPI spec. Study registration/protocol is not completed evidence; sponsor text not peer review. Record registry result publication dates and missing outcome data. | https://clinicaltrials.gov/data-about-studies/learn-about-api ; https://clinicaltrials.gov/api/v2/studies ; https://clinicaltrials.gov/api/oas/v2/ctg-oas-v2.yaml |


### 20.2 Eighteen exact-strain and preparation records

| ID | Exact strain/material and form | Target / setting / observed result | Delivery and commercial identity | Primary source / boundaries |
|---|---|---|---|---|
| TP01 | Bacillus subtilis MB40; OPTI-BIOME MB40; viable spores. Category: pathogen antagonism. | Staphylococcus aureus intestinal and nasal colonization CFU. **human randomized placebo-controlled phase 2 trial**. positive; stool S. aureus reduced 96.8%, nasal 65.4% in trial. | oral spore capsules used; target outcome measured after oral administration; do not infer all B. subtilis spores equivalent. Product: OPTI-BIOME ingredient identity verified in trial; retail SKU not verified. | [PMID:36646104; DOI:10.1016/S2666-5247(22)00322-6](https://pmc.ncbi.nlm.nih.gov/articles/PMC9932624/). decolonization endpoint, not prevention/cure of clinical infection; strain chosen for higher fengycin output. |
| TP02 | Bacillus subtilis 6D1; cell-free extract/lipopeptide fraction. Category: biofilm. | Staphylococcus aureus ATCC 29213 biofilm biomass, mature biofilm disruption, antibiotic sensitization. **in vitro biofilm assays; epithelial cell experiments**. positive; cell-free extract suppressed formation and disrupted mature biofilm; Agr quorum-sensing mechanism investigated. | not_tested for oral gut delivery of active extract; no human efficacy claim. Product: not_verified. | [DOI:10.1128/msystems.00712-24](https://journals.asm.org/doi/10.1128/msystems.00712-24). 6D1 extract not interchangeable with MB40 or generic surfactin; concentration/matrix dependent. |
| TP03 | Limosilactobacillus reuteri DSM 17938 + glycerol; paper: Lactobacillus reuteri 17938; precursor-dependent reuterin production. Category: pathogen antagonism. | Clostridioides difficile CD2015 ribotype 027 growth/invasion. **in vitro and antibiotic-treated human fecal minibioreactor community**. positive conditional on glycerol substrate and reuterin pathway. | community reactor exposure bypasses stomach; no oral human anti-CDI efficacy established by this study. Product: not_verified for combined glycerol formulation. | [PMID:28760934; DOI:10.1128/IAI.00303-17](https://pmc.ncbi.nlm.nih.gov/articles/PMC5607411/). retain glycerol cointervention; do not report comparable-to-vancomycin culture inhibition as comparable clinical treatment. |
| TP04 | Escherichia coli Nissle 1917; EcN; native microcin-producing strain. Category: pathogen/pathobiont competition. | competing Enterobacteriaceae including Salmonella under intestinal inflammation. **mouse inflamed-gut model with genetic microcin controls**. positive; microcins supported competitive suppression during inflammation. | oral mouse model; human gut efficacy for these target organisms untested in this paper. Product: not_verified current SKU. | [PMID:27798599; DOI:10.1038/nature20557](https://www.nature.com/articles/nature20557). inflammation/iron ecology conditions matter; retain target strains from full text before species-specific ranking. |
| TP05 | Bifidobacterium bifidum PRL2010 + Bifidobacterium breve UCC2003; mucin-based coculture. Category: cross-feeding. | growth of B. breve UCC2003 using products released by PRL2010. **in vitro coculture and molecular analysis**. positive cross-feeding. | not_tested oral delivery of paired formulation. Product: not_verified. | [DOI:10.1186/s12866-014-0282-7](https://link.springer.com/article/10.1186/s12866-014-0282-7). ecological support not demonstrated host benefit; mucin substrate should not be relabeled as a generic prebiotic recommendation. |
| TP06 | Bifidobacterium adolescentis 22L or Bifidobacterium breve 12L or Bifidobacterium thermophilum JCM1207 with Bifidobacterium bifidum PRL2010; separate pairwise cocultures, not one four-strain blend. Category: cross-feeding. | PRL2010 growth from starch/xylan degradation products. **in vitro pairwise coculture**. positive with substrate/partner specificity. | not_tested for oral blend. Product: not_verified. | [PMID:26441950; DOI:10.3389/fmicb.2015.01030](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2015.01030/full). import as separate atomic pair/substrate edges after full-text extraction; not proof an arbitrary bifidobacterial mix works. |
| TP07 | Lacticaseibacillus rhamnosus GG-derived p40; paper: Lactobacillus rhamnosus GG; soluble protein, not whole probiotic. Category: barrier mechanism. | ADAM17/HB-EGF/EGFR activation; epithelial survival/barrier pathway. **T84 human cell line, mouse colon epithelial cells, mouse experiments**. positive mechanistic pathway evidence. | not_tested for ordinary oral p40; purified-protein exposure not equivalent to viable LGG supplementation. Product: not_verified for p40 product. | [PMID:24043629](https://pmc.ncbi.nlm.nih.gov/articles/PMC3798544/). protein/activity material identity must remain separate from LGG cells. |
| TP08 | Lactiplantibacillus plantarum WCFS1; paper: Lactobacillus plantarum WCFS1. Category: barrier mechanism. | duodenal ZO-1/occludin localization; chemically stressed epithelial permeability. **human mechanistic intraduodenal experiment plus in vitro epithelium**. positive tight-junction localization and in vitro protection. | direct_site_administration: feeding catheter into duodenum bypassed stomach; does not validate ordinary capsule delivery. Product: not_verified. | [PMID:20224007; DOI:10.1152/ajpgi.00327.2009](https://pubmed.ncbi.nlm.nih.gov/20224007/). healthy-subject tissue biomarker, not cure of a clinical leaky-gut syndrome. |
| TP09 | Bifidobacterium longum 35624; paper calls it Bifidobacterium infantis 35624; encapsulated viable cells. Category: IBS. | abdominal pain and composite symptoms in women with IBS. **human dose-ranging randomized placebo-controlled trial**. positive at tested 10^8 CFU arm; other tested doses not significant. | encapsulated oral outcome observed; largest-dose formulation problems reported; independently verify current final dosage form. Product: not_verified current SKU. | [PMID:16863564; DOI:10.1111/j.1572-0241.2006.00734.x](https://pubmed.ncbi.nlm.nih.gov/16863564/). keep efficacy and formulation/dose specific; do not generalize to B. infantis species. |
| TP10 | Lactiplantibacillus plantarum 299v; Lactobacillus plantarum 299v; DSM 9843; LP299V. Category: IBS. | abdominal pain/bloating and global symptoms. **human randomized placebo-controlled trial**. positive in 2012 trial. | independent human GI survival evidence: PMID 15702859; viable stool recovery supports transit, not permanent colonization. Product: not_verified current SKU. | [PMID:22912552; DOI:10.3748/wjg.v18.i30.4012](https://pubmed.ncbi.nlm.nih.gov/22912552/). one positive trial is not comprehensive strain evidence summary; clinical subtype transfer requires care. |
| TP11 | Bifidobacterium animalis subsp. lactis HN019; paper: Bifidobacterium lactis HN019; viable powder in maltodextrin. Category: functional constipation. | complete spontaneous bowel movements per week. **human multicenter triple-blind randomized trial, 229 participants, 8 weeks**. null primary endpoint; adjusted difference 0.14 CSBM/week (95% CI -0.17 to 0.45), P=.37. | strain DNA recovered/increased in stool is a separate exposure measure, not proof of viable persistence or efficacy. Product: not_verified current SKU. | [PMID:39356506; DOI:10.1001/jamanetworkopen.2024.36888](https://jamanetwork.com/journals/jamanetworkopen/fullarticle/2824333). excellent null-evidence regression test; do not omit because older guide entries suggest constipation benefit. |
| TP12 | Lactobacillus acidophilus CL1285 + Lacticaseibacillus casei LBC80R; paper: L. casei LBC80R; proprietary two-strain formulation. Category: AAD/CDAD prevention. | antibiotic-associated diarrhea and C. difficile-associated diarrhea incidence. **human randomized double-blind dose-response placebo trial**. positive prevention results reported. | oral proprietary formulation outcome; no independent strain survival measurement extracted. Product: Bio-K+ historical studied two-strain formula; current additional-strain products require distinct formulation ID. | [PMID:20145608](https://pubmed.ncbi.nlm.nih.gov/20145608/). mixture-level evidence; cannot allocate effect to either component or assume current three-strain product identical. |
| TP13 | Limosilactobacillus reuteri DSM 17938 + ATCC PTA 6475; paper: Lactobacillus reuteri; two-strain preparation. Category: H. pylori adjunct. | urea breath-test signal and antibiotic-associated adverse events. **human randomized double-blind placebo-controlled trial**. positive bacterial-load/side-effect results; do not infer confirmed stand-alone eradication. | oral studied combination; therapeutic target is gastric so enteric-only delivery could change relevance. Product: BioGaia study page verifies the exact pair; current product/dose equivalence needs label check. | [PMID:24296423; DOI:10.1097/MCG.0000000000000007](https://pubmed.ncbi.nlm.nih.gov/24296423/). adjunct, bacterial load and eradication are separate endpoints. |
| TP14 | Limosilactobacillus reuteri NCIMB 30242; paper: Lactobacillus reuteri NCIMB 30242; bile-salt-hydrolase-active capsule. Category: metabolic. | LDL cholesterol/sterol absorption. **human randomized placebo-controlled trial in hypercholesterolemic adults**. positive lipid endpoint. | oral capsule intervention; separate microencapsulated-yogurt RCT PMID 22067612; delayed-vs-standard release mechanistic pilot PMID 25612224. Product: not_verified current SKU. | [PMID:22990854; DOI:10.1038/ejcn.2012.126](https://pubmed.ncbi.nlm.nih.gov/22990854/). lipid biomarkers not cardiovascular events; delivery formulations should be separate nodes. |
| TP15 | Lacticaseibacillus rhamnosus GG; Lactobacillus rhamnosus GG; ATCC 53103. Category: acute gastroenteritis. | moderate-to-severe gastroenteritis outcomes in preschool children. **human multicenter randomized placebo-controlled trial**. no better outcomes than placebo. | human oral intervention; independent viable colonic-mucosal recovery study PMID 9872808 does not reverse null clinical result. Product: not_verified current SKU. | [PMID:30462938; DOI:10.1056/NEJMoa1802598](https://pubmed.ncbi.nlm.nih.gov/30462938/). maintain clinical setting/population; null for this endpoint does not erase LGG mechanistic or other-indication evidence. |
| TP16 | Bifidobacterium bifidum HI-MIMBb75; SYN-HI-001; heat-inactivated MIMBb75. Category: IBS. | combined abdominal pain and global symptom relief. **human multicenter randomized double-blind placebo-controlled trial**. positive. | viability not required; not_applicable_nonviable for bacterial survival; delivery of intact inactivated material still relevant. Product: not_verified current SKU. | [PMID:32277872](https://pubmed.ncbi.nlm.nih.gov/32277872/). separate from live MIMBb75; do not use live-CFU or engraftment requirement to exclude nonviable candidates. |
| TP17 | Lactiplantibacillus plantarum 299v; DSM 9843; viable preparation. Category: delivery. | viable gastrointestinal transit/fecal recovery with and without acid inhibition. **human placebo-controlled double-blind study**. positive survival through GI transit. | viable_fecal_recovery; supports delivery plausibility, not exact site concentration or durable engraftment. Product: not_verified current SKU. | [PMID:15702859](https://pubmed.ncbi.nlm.nih.gov/15702859/). link to TP10 as separate evidence; formulation equivalence needs verification. |
| TP18 | Limosilactobacillus reuteri NCIMB 30242; delayed-release vs standard-release capsule. Category: delivery/formulation. | bile acids, FGF-19, sterol absorption. **human small mechanistic randomized pilot, 10 adults**. formulation-dependent biochemical response. | engineered delayed-release formulation tested; mechanistic exposure evidence, not general proof enteric coating improves all strain outcomes. Product: not_verified current SKU. | [PMID:25612224](https://pubmed.ncbi.nlm.nih.gov/25612224/). small pilot with surrogate markers; capsule design distinct from strain identity. |

**Required strain-specific counterevidence for TP04:** Nissle 1917 has published colibactin-dependent epithelial DNA cross-linking/mutagenicity experiments ([2021 primary study, PMID 34378987](https://pubmed.ncbi.nlm.nih.gov/34378987/)), alongside an earlier study reporting no detected genotoxicity under its different assays ([2020, PMID 32363034](https://pubmed.ncbi.nlm.nih.gov/32363034/)). Preserve assay/context differences and the favorable microcin observations together. Do not label a probiotic strain universally safe or harmful from its species name, or turn these experiments into a numerical human cancer risk. A [2019 mechanistic study](https://doi.org/10.1371/journal.ppat.1008029) separates ClbP peptidase and antibacterial activity, but engineered derivatives are different strain entities. A [2021 oligosaccharide experiment](https://pubmed.ncbi.nlm.nih.gov/33596864/) also links particular tested prebiotics to colibactin expression/genotoxicity; extract exact preparations before binding feeding tradeoffs. These are additional observation edges under TP04, not independent positive efficacy counts.

### 20.3 Additional primary-source extraction requirements

- **Saccharomyces boulardii CNCM I-745 → C. difficile biofilm formation**, primary in vitro study DOI https://doi.org/10.3390/microorganisms10061082 (2022), https://www.mdpi.com/2076-2607/10/6/1082 . This is explicitly the I-745 strain. Extract direction, target isolates, biofilm assay and preformed-versus-formation distinction before adding a quantitative edge.
- **S. boulardii CNCM I-745 → epidemic C. difficile strain challenge**, primary 2016 study https://pmc.ncbi.nlm.nih.gov/articles/PMC5142203/ . Distinguish strain-labeled primary study from older S. boulardii protease articles whose commercial/strain identity needs source verification. Older human-colonic-mucosa ex vivo protease study PMID 9864230 https://pubmed.ncbi.nlm.nih.gov/9864230/ is an additional mechanistic source; do not attach I-745 solely from species name.
- **L. plantarum WCFS1, CIP104448, TIFN101 → NSAID-stressed small-intestinal barrier**, separate oral human crossover periods: PMID 28045137 https://pubmed.ncbi.nlm.nih.gov/28045137/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC5206730/ . This oral study should be extracted alongside TP08 because transcriptional changes and functional permeability benefit are not interchangeable.
- **B. lactis HN019 2025 8-week functional constipation trial**, DOI https://doi.org/10.1002/mnfr.70081 . Recent study must be read and deduplicated against the 2024 trial before declaring the whole HN019 evidence base positive or negative.
- **L. reuteri DSM 17938 + ATCC PTA 6475 + PPI without antibiotics**, 2019 RCT PMID 31065301 https://pubmed.ncbi.nlm.nih.gov/31065301/ ; full https://pmc.ncbi.nlm.nih.gov/articles/PMC6466906/ . Useful explicit limit on extrapolating the adjunctive H. pylori signal to stand-alone eradication.
- **MIMBb75 live cells**, 2011 IBS RCT PMID 21418261 https://pubmed.ncbi.nlm.nih.gov/21418261/ ; separate node from TP16 heat-inactivated formulation.

### 20.4 Integration requirements

1. Import Probio-Ichnos JSON as a reversible source snapshot retaining PMID/property flags and license. Enrich antimicrobial/biofilm targets from full texts instead of returning a generic antimicrobial boolean.
2. Join BacDive/NCBI identifiers with explicit alias evidence; retain verbatim name so taxonomic renaming cannot drop studies. Strain similarity is not automatic therapeutic equivalence.
3. Ingest IPDB as predicted capabilities, not an efficacy score. Positive mechanistic, null human, and unmeasured delivery edges can coexist for the same strain.
4. Attach human-guide references and registry outcomes after broad discovery, not as eligibility gates. UK guide status is currently a freshness warning, not a live product authority.
5. Model retail product, historical trial formulation, viable strain, inactivated cells, purified molecule, and cell-free extract as distinct entities. Equivalence edges require sourced identity, amounts, viability, release/matrix and date.
6. Return broad target-matched options with explicit evidence tier, tested material, reachable compartment, and main missing link. User-requested all-evidence retrieval is compatible with distinguishing hypothesis from established benefit.


## 21. Barrier, HMO, prebiotic and ecological-tradeoff records

Records BP01–BP25, their product-identity checks and opposing observations extend the same intervention graph. They must not create a second scoring engine.

### 21.1 Decisions the recommender should make explicit

- Include zinc carnosine/polaprezinc, glutamine, BPC-157, KPV, individual HMO structures, and recombinant human lactoferrin as distinct candidates. Evidence level determines confidence and presentation, not whether a mechanistically credible candidate can be represented.
- Represent `clinical benefit`, `barrier function`, `inflammation`, `taxon abundance`, `substrate use`, `cross-feeding`, `virulence`, `tolerability`, `product identity`, and `long-term safety` separately. One cannot substitute for another.
- Stool DNA is a compositional and functional-potential measurement. It does not directly measure epithelial permeability, mucus thickness, host tight-junction integrity, bacterial viability, or metabolite flux. Keep microbial barrier hypotheses separate from directly measured barrier dysfunction. The cited permeability trials used orally administered sugar probes; the experimental HMO work used tracer permeability/cell models.
- Microbiome data at species level cannot establish a particular strain's HMO utilization, mucolysis, inflammatory polysaccharide, or toxin phenotype. A substrate-use capability does not establish that oral supplementation increases the organism in a complex human community, or that the net clinical result is harmful.
- Explicitly model opposing paths. A substrate can support a desired organism, a cross-fed organism, and host repair simultaneously. The net outcome depends on structure, strain, concentration, community partners, inflammation, site, host absorption, and formulation.
- There is no retrieved direct evidence that oral glutamine selectively causes an R. gnavus bloom. Do not turn generic nitrogen use, glutamine-dependent genes, or host glutamine association into a glutamine contraindication. There is direct context-dependent L-serine/AIEC evidence, but that is not glutamine/R. gnavus evidence.

### 21.2 Twenty-five concrete edges

| ID | Intervention → target / observed direction | Evidence and context | Tradeoff, resolution, or boundary | Primary source |
|---|---|---|---|---|
| BP01 | Oral zinc carnosine → less NSAID-induced small-intestinal hyperpermeability | Randomized crossover, 10 healthy volunteers; indomethacin increased lactulose:rhamnose permeability on placebo, not during zinc carnosine coadministration. | A small drug-challenge trial supports this injury context. It does not validate stool-DNA diagnosis or all causes of intestinal symptoms. | [Mahmood 2007, PMID 16777920](https://pubmed.ncbi.nlm.nih.gov/16777920/), DOI 10.1136/gut.2006.099929. |
| BP02 | Zinc carnosine, alone/with bovine colostrum → attenuated exercise-associated permeability | Four-arm double-blind crossover in 8 volunteers; exercise challenge; parallel epithelial-cell experiments. | Small, short study. Combination effects and heat/exercise context should not be transferred wholesale to chronic disease or zinc alone. | [Davison 2016, PMID 27357095](https://pubmed.ncbi.nlm.nih.gov/27357095/), DOI 10.3945/ajcn.116.134403. |
| BP03 | Rectal polaprezinc → mucosal-healing signal in active ulcerative colitis | Investigator-blinded randomized study, 28 patients; adjunctive enema, short observation. | Preserve route and add-on setting: this is not evidence that an oral supplement reproduces local enema exposure. | [Itagaki 2014, PMID 24286534](https://pubmed.ncbi.nlm.nih.gov/24286534/). |
| BP04 | Oral glutamine → improved IBS symptoms and sugar-probe permeability | Double-blind placebo-controlled trial in postinfectious IBS-D selected for increased permeability; 106 completers over 8 weeks. Primary symptom response: 79.6% versus 5.8%. | Stronger match for this particular phenotype; replication/generalization remain limitations. It does not show glutamine is necessary for every dysbiosis profile. | [Zhou 2019, PMID 30108163](https://pubmed.ncbi.nlm.nih.gov/30108163/). |
| BP05 | Glutamine → no demonstrated superiority for permeability in Crohn's disease | Small randomized adult study; no significant permeability benefit. A separate pediatric active-Crohn trial found no advantage for glutamine-enriched versus standard polymeric diet. | Negative condition-specific evidence limits a universal barrier-healing label. Do not collapse Crohn's and postinfectious IBS-D. | [Den Hond 1999, PMID 9888411](https://pubmed.ncbi.nlm.nih.gov/9888411/); [Akobeng 2000, PMID 10630444](https://pubmed.ncbi.nlm.nih.gov/10630444/). |
| BP06 | Glutamine → reduced EHEC type-III-secretion virulence and intestinal colonization in models | Bacterial experiments plus mouse infection: nitrogen-metabolic signaling repressed T3SS; host defense also improved. | Useful counterexample to “an amino acid can be metabolized, therefore it worsens pathogens.” Not a proven human EHEC treatment. | [Microbiology Spectrum 2023](https://journals.asm.org/doi/10.1128/spectrum.00975-23), DOI 10.1128/spectrum.00975-23. |
| BP07 | L-glutamine + suitable germinant environment → C. difficile spore cogermination | In-vitro strain-specific cogerminant study; glutamine was among recognized amino acids in strain UK1. | A real possible microbial interaction, but germination is not equivalent to outgrowth, toxin disease, or net effect of orally ingested glutamine. Keep bile-acid context and tested strain. | [Hierarchical recognition of amino acid co-germinants, 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC5844826/). |
| BP08 | Alanyl-glutamine → epithelial restitution after TcdB and improved mouse CDI recovery | Cell and mouse experiments using the dipeptide formulation. | Host repair may oppose the preceding nutrient-use pathway; formulation differs from free glutamine. Human efficacy should not be asserted from this experiment or a trial protocol. | [Intestinal epithelial restitution after TcdB challenge, PMID 23359592](https://pubmed.ncbi.nlm.nih.gov/23359592/). |
| BP09 | BPC-157 → reduced experimental cysteamine colitis injury | Primary rat studies of induced colon lesions, including delayed treatment in chronic injury models. | Eligible as preclinical candidate. Chemical injury is not human IBD; human gut efficacy has not been established by these papers. | [PMID 11595448](https://pubmed.ncbi.nlm.nih.gov/11595448/), [PMID 11595451](https://pubmed.ncbi.nlm.nih.gov/11595451/). |
| BP10 | BPC-157 → improved intestinal anastomosis/fistula healing in rats | Primary ileoileal anastomosis study and colitis-complicated colon-anastomosis study. | Surgical tissue repair is a separate endpoint from human permeability. Article titles referring to clinical development are not themselves human efficacy evidence. | [PMID 17713731](https://pubmed.ncbi.nlm.nih.gov/17713731/), [Klicek 2013 full paper](https://jpp.krakow.pl/journal/archive/10_13/articles/09_article.html). |
| BP11 | KPV uptake through PepT1 → reduced NF-κB/MAPK activation and IL-8 | Human epithelial/T-cell lines; transporter competition and expression experiments support PepT1 dependence. | Mechanistic rationale is more specific than “anti-inflammatory peptide.” Cell exposure is not established oral human exposure; competing peptides mattered experimentally. | [Dalmasso 2008, PMID 18061177](https://pmc.ncbi.nlm.nih.gov/articles/PMC2431115/), DOI 10.1053/j.gastro.2007.10.026. |
| BP12 | Oral KPV → reduced DSS/TNBS colitis inflammation in mice | Drinking-water intervention reduced inflammatory and histological measures in two injury models. | Include as animal efficacy evidence; it does not establish clinical efficacy, a human dosing regimen, or long-term safety. | [Dalmasso 2008](https://pmc.ncbi.nlm.nih.gov/articles/PMC2431115/). |
| BP13 | 2′FL and/or LNnT → increased fecal Bifidobacterium in healthy adults | Placebo-controlled randomized study, 100 adults, 2 weeks; mixtures and individual structures tested. | Primarily microbiome/tolerability endpoints. Response was not universal; increased relative abundance does not establish symptomatic benefit or exclusion of all other utilizers. | [Elison 2016, PMID 27719686](https://pmc.ncbi.nlm.nih.gov/articles/PMC5082288/), DOI 10.1017/S0007114516003354. |
| BP14 | Defined 4:1 2′FL:LNnT → increased fecal Bifidobacterium in IBS without overall symptom aggravation | Randomized placebo-controlled trial, 61 randomized/58 completers, 4 weeks; the higher tested arm increased bifidobacteria. | This study primarily establishes microbiome change/tolerability, not robust placebo-controlled clinical symptom efficacy. Later sample analysis is the same cohort, not independent replication. | [Iribarren 2020, PMID 32536023](https://pubmed.ncbi.nlm.nih.gov/32536023/); [2021 follow-up, PMID 34836092](https://pmc.ncbi.nlm.nih.gov/articles/PMC8622683/). |
| BP15 | 2′FL → bifidobacterial bloom in older adults; immune endpoint null | Six-week randomized study of 89 older adults; primary cytokine-response endpoint was not met, while microbiota and some secondary outcomes changed. | A contemporary example where target engagement is real but the primary desired host effect is not demonstrated. Do not convert secondary responder analyses into a general cognitive/immune benefit. | [Carter 2025, PMID 40738103](https://pmc.ncbi.nlm.nih.gov/articles/PMC12432366/). |
| BP16 | Fermented 2′FL / 2′FL+LNnT products → lower epithelial tracer permeability | SHIME adult fecal fermentation, Caco-2 monolayers, human biopsy-derived gut-on-chip models; barrier-related gene expression and inflammatory changes. | This is fermentation-supernatant/ex-vivo evidence, not a human oral permeability trial. The donor community and fermentation products mediate the effect. | [Šuligoj 2020, PMID 32933181](https://pmc.ncbi.nlm.nih.gov/articles/PMC7551690/), DOI 10.3390/nu12092808. |
| BP17 | HMO structure → strain-specific R. gnavus growth | Both E1 and ATCC29149 grew on 2′FL and 3FL; neither grew on LNT or LNnT; ATCC29149, unlike E1, used 3′SL. | Directly defeats “all HMOs only feed good bacteria.” Equally, no growth in two isolates does not prove universal LNnT immunity to utilization/cross-feeding. No human adverse outcome established. | [Crost 2013](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0076341), DOI 10.1371/journal.pone.0076341. |
| BP18 | Extracellular 2′FL cleavage by community members → enables B. breve utilization | Infant microbiome cultivation/metagenomics found B. breve could use 2′FL in a consortium; metagenomics implicated extracellular fucosidases of coexisting members such as R. gnavus. | Marks beneficial cross-feeding and imperfect exclusivity. The implicated R. gnavus contribution is a metagenomic mechanistic inference, not proven adult clinical treatment response. | [PMID 37973815](https://pubmed.ncbi.nlm.nih.gov/37973815/). |
| BP19 | Resistant potato starch → R. bromii/other degraders → higher fecal butyrate in responders | Human dietary studies show heterogeneous responses; degradation partners and E. rectale abundance predicted response. Potato, maize starch, and inulin were not interchangeable. | R. bromii is not R. gnavus. Fecal butyrate is an endpoint of production minus utilization/absorption, not direct mucosal exposure. Taxonomic response alone does not guarantee functional response. | [Baxter 2019, PMID 30696735](https://pmc.ncbi.nlm.nih.gov/articles/PMC6355990/); [Venkataraman 2016, PMID 27357127](https://pmc.ncbi.nlm.nih.gov/articles/PMC4928258/). |
| BP20 | R. bromii starch breakdown → indirect R. gnavus feeding | Co-culture experiments: R. gnavus ATCC29149 did not use potato soluble/maize resistant starch alone but used products liberated by R. bromii L2-63. | Precisely why a monoculture “cannot use” result does not mean a food cannot feed that organism in a community. Net human harm is untested, and desired butyrate cross-feeding may coexist. | [Crost 2018, PMID 30455672](https://www.frontiersin.org/journals/microbiology/articles/10.3389/fmicb.2018.02558/full), DOI 10.3389/fmicb.2018.02558. |
| BP21 | Inulin → higher R. gnavus representation in a defined community | July 2026 six-species synthetic-community study: inulin added to BHIS medium promoted R. gnavus and reduced several other community members. | Context is a simplified, nutrient-rich batch-culture model. This is a conditional tradeoff signal, not evidence that inulin always increases R. gnavus or causes clinical harm. | [Yu 2026 mSphere](https://journals.asm.org/doi/10.1128/msphere.00289-26), DOI 10.1128/msphere.00289-26. |
| BP22 | Fiber deprivation → increased host-mucus consumption → mucus erosion/pathogen susceptibility | Gnotobiotic mice with a defined human gut community; fiber deprivation shifted community substrate use and increased C. rodentium susceptibility. | Cautions against indiscriminate starvation of the community as a way to suppress a taxon. It does not imply every isolated fiber repairs mucus in all people. | [Desai 2016, PMID 27863247](https://pmc.ncbi.nlm.nih.gov/articles/PMC5131798/), DOI 10.1016/j.cell.2016.10.043. |
| BP23 | Dietary L-serine → competitive advantage for AIEC LF82 in the inflamed gut | Mechanistic bacterial and mouse work: inflammation changed metabolic programs; removing serine blunted blooms in the tested setting. | Keep strain, inflammatory context, and nutrient separate. This does not establish that glutamine feeds R. gnavus or justify general amino-acid/protein restriction. | [Kitamoto 2019 Nature Microbiology, PMID 31686025](https://www.nature.com/articles/s41564-019-0591-6). |
| BP24 | L-serine restriction + mucolytic community → AIEC access to host nutrients/epithelial niche | Follow-up mouse work found A. muciniphila-mediated mucus degradation permitted AIEC to obtain host-derived nutrients during serine restriction; restriction could increase E. coli in that setting. | A documented direction reversal from community context. The right model is conditional rules, not permanent “serine bad” or “Akkermansia good” scores. | [Cell Reports 2022](https://pmc.ncbi.nlm.nih.gov/articles/PMC10903618/), article S2211124722008956. |
| BP25 | Oral recombinant human lactoferrin → attenuated NSAID-induced small-intestinal permeability | Double-blind placebo-controlled crossover in 15 healthy volunteers using an indomethacin challenge. | Human barrier evidence exists for this ingredient class, but the 2003 recombinant product is not automatically equivalent to effera/Kepos in manufacturing, glycosylation, quantity, or finished formulation. | [Troost 2003, PMID 14647223](https://pubmed.ncbi.nlm.nih.gov/14647223/), DOI 10.1038/sj.ejcn.1601727. |

### 21.3 Kepos / “Keppos” resolution

Likely intended product: **kēpos for Digestion**, [official product page](https://trykepos.com/products/kepos-for-digestion-daily-prebiotic), opened 2026-09-22 . Its listed ingredients are **kpHMO™ (human milk oligosaccharides)** and **effera® (human milk lactoferrin)**. It is fermentation-produced, not collected human milk. These are product identity claims, not proof of efficacy.

The [official buyer-guide FAQ](https://trykepos.com/blogs/blogs/best-hmo-supplements-how-to-choose-the-right-one-for-you) explicitly states that the specific HMO structures are not disclosed (, paragraph around line 133). It claims neutral, fucosylated and sialylated category coverage. Earlier brand pages mention 2′FL, but that is insufficient to lock the current product to pure 2′FL or to the 4:1 2′FL/LNnT clinical mixture.

**Product mapping result:** manufacturer and two ingredient families verified; individual HMO identities/ratios and equivalence to tested formulations unresolved. Do not fill these gaps with a guess. Ingredient-family evidence can remain visible, with the branded product's match confidence reduced and a label/lot composition requirement before structure-specific ranking. The marketed claim that harmful microbes cannot access HMOs is contradicted as a universal substrate-utilization statement by BP17/BP18. That contradiction is not proof that the finished product is harmful.

**Current effera primary data:** [Peterson et al., 2026, PMID 42178844](https://pubmed.ncbi.nlm.nih.gov/42178844/)  reports exploratory microbiome outcomes from a randomized double-blind study of 66 healthy adults, two effera arms versus bovine lactoferrin for 28 days. Diversity had no main treatment effect; genus/metabolite changes occurred, including Faecalibacterium increase in the higher effera arm. These were explicitly hypothesis-generating; there was no placebo arm. This is not a finished-Kepos efficacy trial or a demonstrated gut-permeability benefit. The [effera barrier trial NCT07035964](https://clinicaltrials.gov/study/NCT07035964) exists, but usable results/status were not exposed by the page retrieval, so none are asserted here.

### 21.4 Separate safety, exposure, and quality axes

**Zinc carnosine:** Count elemental zinc from all products, not chelate mass. Chronic excess can impair copper absorption. [NIH ODS zinc fact sheet](https://ods.od.nih.gov/factsheets/Zinc-HealthProfessional/)  gives an adult tolerable upper intake of 40 mg/day outside clinician-directed treatment and discusses copper risk with sustained high intake. This is a known dose/exposure risk, not an argument against representing the intervention.

**Glutamine:** Do not say harmless in every setting. The [REDOXS randomized trial, PMID 23594003](https://pubmed.ncbi.nlm.nih.gov/23594003/) included 1,223 ventilated ICU patients with multiorgan failure; early supplementation was associated with increased mortality. This is a materially different high-risk population and exposure context from the IBS trial and should trigger a separate restriction, not contaminate inference about usual outpatient glutamine or microbial feeding. Acute clinical context outranks a stool abundance heuristic.

**BPC-157:** [Xu 2020 preclinical toxicology, PMID 32334036](https://pubmed.ncbi.nlm.nih.gov/32334036/) found no serious toxicity in studied mice, rats, rabbits and dogs and reported no genetic or embryo-fetal toxicity. [He 2022 PK, PMID 36588717](https://pubmed.ncbi.nlm.nih.gov/36588717/) characterized IV/IM exposure in rats/dogs; those data do not establish oral human bioavailability. A [2025 IV human pilot, PMID 40131143](https://pubmed.ncbi.nlm.nih.gov/40131143/) reported tolerability in only two adults, with no efficacy conclusion possible. Thus “no human data at all” is inaccurate, while “established safe/effective for gut repair” is also inaccurate. Separate (a) animal efficacy, (b) short-term animal toxicology, (c) tiny human exposure observations, (d) human gut efficacy unknown, and (e) formulation/identity/sterility/impurity risks.

**KPV:** The [FDA 2026 primary evaluation](https://www.fda.gov/media/193346/download)  identified preclinical pharmacology but no human clinical/PK evidence, no formal acute/repeat-dose/genotoxic/reproductive toxicology package, and insufficient immunogenicity/aggregation information. Searches through December 2025 found no adverse-event cases, which cannot establish safety. These are information/quality gaps, not evidence of demonstrated clinical toxicity. Keep the positive cell/mouse edges BP11/BP12. The briefing is an evaluation/proposal, not itself a final legal determination. Do not conflate topical formulation review with proof of oral delivery or compound salts/derivatives with identical clinical exposure.

### 21.5 Additional precision checks and implementation notes

1. Preserve source type even when the paper uses human-derived cells: `human cell line` and `human organoid ex vivo` are not `human trial`.
2. Distinguish `clinical endpoint not measured`, `measured but null`, `evidence not found`, and `harm observed`. These require different recommender behavior.
3. A negative R. gnavus monoculture edge for LNnT should use `no_growth_in_tested_strains_conditions`, not `does_not_feed_pathobionts`.
4. “Ruminococcus” genus-wide actions are especially unsafe: R. bromii and R. gnavus have different substrate roles, and historical taxonomy has changed. Preserve aliases and exact strain identifiers.
5. The [2024 R. gnavus arginine/NO study, PMID 39089585](https://pubmed.ncbi.nlm.nih.gov/39089585/) identifies a host–microbe arginine–NOS2 feedback circuit. It also reports elevated glutamine in spent R. gnavus culture medium. It is not a glutamine-supplement feeding trial; it does not support assigning a harmful oral-glutamine→R. gnavus edge..
6. The commonly cited [Park 2020 BPC-157 NSAID/permeability article, PMID 32445447](https://pubmed.ncbi.nlm.nih.gov/32445447/) describes tight-junction/cytoprotection mechanisms, but its publisher abstract labels it a review. Do not count it as an independent randomized human barrier trial or silently inflate primary-evidence count. Primary animal healing studies are provided in BP09/BP10.
7. Appropriate recommendation card fields: `candidate`, `exact_formulation`, `target/problem`, `mechanistic_path`, `organism/strain`, `community_conditions`, `human_population`, `route`, `endpoint`, `direction`, `evidence_design`, `source`, `effect_uncertainty`, `opposing_edges`, `safety_observed`, `safety_unknown`, `product_equivalence`, and `what_would_change_ranking`.
8. Barrier candidates should rank higher when the clinical phenotype matches the study (e.g. selected postinfectious IBS-D, NSAID injury) and lower when only a stool-taxonomy hypothesis matches. BPC/KPV should remain discoverable as investigational mechanisms with unknown human gut efficacy; neither should receive an invented safe regimen.
9. The [HMO open-label IBS study, PMID 33512807](https://pubmed.ncbi.nlm.nih.gov/33512807/) reports substantial symptom improvement but lacks placebo control. It can support a human signal while remaining below randomized symptom efficacy. It is distinct from the randomized bifidogenicity trial; do not merge their outcome claims.



## 22. Targeted antimicrobial and biofilm intervention records

TM01–TM30 are mandatory seed coverage across distinct endpoint types. TM25 and TP03 describe the same reuterin study: maintain both index aliases but one set of underlying observations, not duplicate replication. Apply this deduplication rule to overlaps with older I/SYN records as well.

### 22.1 Endpoint and delivery labels

- `PLANKTONIC`: MIC/MBC/time-kill or growth inhibition in free-living cells. Not a biofilm result.
- `BIOFILM_PREVENTION`: compound present before/during attachment or development. Cannot inherit mature-biofilm clearance.
- `ESTABLISHED_BIOMASS`: residual crystal-violet (CV), imaging, or matrix after treating a preformed biofilm. CV loss is not bacterial killing.
- `ESTABLISHED_VIABILITY`: CFU, metabolic assay, or membrane staining after treatment; preserve which assay. Metabolic suppression is not necessarily loss of culturability.
- `DISPERSAL`: movement/loss from established attached biomass; viable released cells may remain.
- `ANIMAL_COLONIZATION_PREVENTION`, `HUMAN_CARRIAGE_CLEARANCE`, `HUMAN_ERADICATION_TEST`, `COMMUNITY_ABUNDANCE`: distinct outcomes, never silently converted into one another.
- Unless an oral route is stated below, delivery is **direct exposure in culture**, with human gut-lumen exposure, acid survival, enteric release, mucosal concentration, microbiome selectivity and tolerability **unestablished in that study**. Store these fields as unknown, rather than assuming failure or success. Laboratory concentrations below describe experiments, not proposed doses.

### 22.2 Thirty primary seed edges

Each table row is a candidate edge (some rows deliberately retain two separately measured endpoints for later splitting). Repeated rows sharing a paper do not represent independent replication. `A` = abstract/primary indexed excerpts sufficient for the stated narrow claim; `M` = relevant methods/results were retrieved. Unknown strain IDs or numeric values must remain null until full-text extraction.

| ID | Exact target and preparation | Measured endpoint and result | Model, delivery and limits | Primary identifier / access |
|---|---|---|---|---|
| TM01 | Uropathogenic *E. coli* CFT073/J96 ← synthesized allicin, HPLC purity 96% | `BIOFILM_PREVENTION`: 12–50 µg/mL, growth-normalized CV; 50 µg/mL lowered CFT073 attached biomass about 33%, J96 about 17%. | M; urinary-pathogen laboratory model. Do not relabel as intestinal biofilm. Methods specify synthetic allicin although abstract describes garlic origin. | Yang et al. 2016, [DOI 10.3390/ijms17070979](https://doi.org/10.3390/ijms17070979); [PMID 27367677 / full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC4964365/). |
| TM02 | UPEC CFT073/J96 ← same synthetic allicin | `ESTABLISHED_BIOMASS` / author-described `DISPERSAL`: 50 µg/mL for 1 h lowered preformed attached CV signal about 40%/30%. | M; separate treatment-after-formation experiment. Biomass dispersal does not establish killing or eradication. | Same Yang study as TM01. |
| TM03 | UPEC ← carvacrol-rich oregano essential oil | `BIOFILM_PREVENTION`: inhibition at subinhibitory concentrations below 0.01%; curli/virulence assays also studied. | A; urinary-pathogen surrogate; exact chemotype/lot matters. No intestinal delivery experiment. | Lee et al. 2017, [DOI 10.1111/jam.13602](https://doi.org/10.1111/jam.13602); [PMID 28980415](https://pubmed.ncbi.nlm.nih.gov/28980415/). |
| TM04 | UPEC ← carvacrol, isolated compound | `BIOFILM_PREVENTION`: reduced formation at subinhibitory concentrations below 0.01%. | A; component tested in same screen as TM03. Store separately from oregano oil. | Same Lee study as TM03. |
| TM05 | UPEC ← thymol, isolated compound | `BIOFILM_PREVENTION`: reduced formation at subinhibitory concentrations below 0.01%. | A; same study also tested thymol-rich thyme-red oil; that is a separate preparation, not proof all thyme products work. | Same Lee study as TM03. |
| TM06 | *K. pneumoniae*, strong-biofilm clinical isolates ← berberine | `BIOFILM_PREVENTION`: reported MBIC 0.0635 mg/mL; seven strong producers selected from 35 isolates. | A; plate-based isolates, no intestinal colonization endpoint. Preserve MBIC definition when extracting full methods. | Magesh et al. 2013, [PMID 24377137](https://pubmed.ncbi.nlm.nih.gov/24377137/), *Indian J Exp Biol* 51:764–772. No DOI confirmed. |
| TM07 | Carbapenem-resistant *K. pneumoniae* ← eugenol | `PLANKTONIC`, `BIOFILM_PREVENTION`, and biofilm-associated cell inactivation reported; inhibition of biofilm-associated gene expression. | A; compound rather than clove oil. Extract mature-biofilm assay, concentrations and CFU/metabolic readout before assigning a quantitative mature-biofilm effect. No gut model. | Qian et al. 2020, [DOI 10.1016/j.micpath.2019.103924](https://doi.org/10.1016/j.micpath.2019.103924); [PMID 31837416](https://pubmed.ncbi.nlm.nih.gov/31837416/). |
| TM08 | *E. faecalis*, ten biofilm-forming clinical isolates from UTI patients ← berberine hydrochloride | `BIOFILM_PREVENTION`: reduced biofilm development, alongside reduced sortase A/esp expression. | A; strains selected from 99 urine isolates. Urinary source, not evidence for intestinal mucosal biofilms. Exact tested concentrations unresolved from abstract. | Chen et al. 2016, [DOI 10.1016/j.micres.2016.03.003](https://doi.org/10.1016/j.micres.2016.03.003); [PMID 27242142](https://pubmed.ncbi.nlm.nih.gov/27242142/). |
| TM09 | Same *E. faecalis* isolates ← berberine hydrochloride | `DISPERSAL`: the study separately reports enhanced dispersion of established biofilm. | A; retain author-reported endpoint; do not turn dispersion into complete killing. Quantitative assay extraction pending. | Same Chen study as TM08. |
| TM10 | *C. albicans* clinical isolate 475/15 ← *Artemisia absinthium* herb ethanol extract | `BIOFILM_PREVENTION`: >50% CV reduction at MIC, 0.5 mg/mL, during 24 h development. | M; extracted with absolute ethanol, dried and redissolved in 30% ethanol. Not wormwood tea, essential oil, *A. annua*, or purified artemisinin. | Ivanov et al. 2021, [DOI 10.1155/2021/9961089](https://doi.org/10.1155/2021/9961089); [PMID 34335850 / full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC8324356/). |
| TM11 | *C. krusei* H1/16 ← same *A. absinthium* ethanol extract | `BIOFILM_PREVENTION`: >50% CV reduction at 0.5 MIC; planktonic MIC 1 mg/mL. | M; in-vitro clinical isolate. No ingestion or gut-colonization evidence. | Same Ivanov study as TM10. |
| TM12 | *K. pneumoniae* ATCC 13883 ← same *A. absinthium* extract | `BIOFILM_PREVENTION`: <30% CV reduction even at MIC. | M; weak result belongs in registry. Candida effect cannot be copied onto Klebsiella. | Same Ivanov study as TM10. |
| TM13 | *K. pneumoniae*, experimentally meropenem-resistant strains ← allicin **plus meropenem** | `PLANKTONIC_SYNERGY`: checkerboard, time-kill and isobologram assays; combination improved activity relative to meropenem alone. Also mouse pneumonia benefit. | A; direct Klebsiella evidence. **Not a biofilm study, not gut decolonization, not allicin monotherapy.** Exact allicin preparation/strain identifiers and in-vivo route need full text. | Liu et al. 2026, [DOI 10.1111/apm.70164](https://doi.org/10.1111/apm.70164); [PMID 41703979](https://pubmed.ncbi.nlm.nih.gov/41703979/). |
| TM14 | Carbapenem-resistant UPEC ← thymoquinone | `PLANKTONIC`: antibacterial effect at 256 µg/mL. `BIOFILM_PREVENTION` and post-formation effects also reported with CV assays. | M, partial; urinary isolates. Keep antibacterial concentration separate from biofilm concentrations. Author word “eradication” must not become complete viable-cell eradication without assay/LOD verification. | Jin & Eom 2024, [DOI 10.1007/s12088-024-01231-8](https://doi.org/10.1007/s12088-024-01231-8); [PMID 39678958 / full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC11645355/). |
| TM15 | *C. glabrata* oral isolates from ICU patients ← thymoquinone | `PLANKTONIC`: MIC50 50 µg/mL; `BIOFILM_PREVENTION`: roughly twofold reduction in formation at MIC50, with lower EPA6 expression. | M, partial; oral isolates, not intestinal colonization. Purified compound, not black-seed oil equivalence. | Nouri et al. 2023, [DOI 10.3390/metabo13040580](https://doi.org/10.3390/metabo13040580); [PMID 37110238 / full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC10143056/). |
| TM16 | Two *C. difficile* strains ← *Nigella sativa* black-seed oil | `PLANKTONIC_GROWTH`: agar-disc inhibition, strongest at 2% oil. | M, partial; oil diluted in methanol with solvent control; water extract also tested. **No biofilm, spores, oral delivery, human efficacy or eradication endpoint.** | Aljarallah 2016, [DOI 10.1016/j.jtumed.2016.05.006](https://doi.org/10.1016/j.jtumed.2016.05.006); [primary full-text PDF](https://applications.emro.who.int/imemrf/J_Taibah_Univ_Med_Sci/J_Taibah_Univ_Med_Sci_2016_11_5_427_431.pdf). |
| TM17 | *H. pylori* infection in adults with non-ulcer dyspepsia ← ground *N. sativa* seed capsules **plus omeprazole** | `HUMAN_ERADICATION_TEST`: stool-antigen-negative 4 weeks after treatment; 2 g/day seed arm 66.7%, triple-therapy comparator 82.6%; 1 g and 3 g arms lower. | M; 88 participants, open trial; methods use alternate-subject allocation. Four weeks oral seed, PPI cointervention. Lack of significant difference is **not demonstrated equivalence/noninferiority**. No biofilm assay. | Salem et al. 2010, [PMID 20616418 / full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC3003218/). |
| TM18 | Refractory *H. pylori* infection ← oral NAC pretreatment followed by culture-guided antibiotics | `HUMAN_ERADICATION_TEST`: 13/20 versus 4/20, after ≥4 previous failures; gastric biofilm assessment and in-vitro biofilm work also performed. | M, partial; open randomized trial, NAC 600 mg/day for 7 days before antibiotics. Combination strategy, not NAC monotherapy; stomach, not colonic biofilm. | Cammarota et al. 2010, [DOI 10.1016/j.cgh.2010.05.006](https://doi.org/10.1016/j.cgh.2010.05.006); [PMID 20478402](https://pubmed.ncbi.nlm.nih.gov/20478402/). |
| TM19 | Treatment-naïve *H. pylori* infection ← NAC concurrently with clarithromycin/amoxicillin/PPI triple therapy | `HUMAN_ERADICATION_TEST`, **null superiority**: adding NAC did not improve first-line eradication. | M, partial; 680 participants, multicenter open randomized trial; NAC 600 mg twice daily, 14 days. Different population/timing from TM18. Preserve both findings. | Chen et al. 2020, [DOI 10.1177/1756284820927306](https://doi.org/10.1177/1756284820927306); [PMID 32821287 / full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC7406927/). |
| TM20 | Mature multispecies biofilm containing *E. faecalis*, *A. naeslundii*, *L. salivarius*, *S. mutans* ← NAC | `ESTABLISHED_BIOMASS` and `ESTABLISHED_VIABILITY`: removal and killing of 3-week biofilms on human dentin; NAC 25, 50, 100 mg/mL tested. | M; intracanal medicament model. Community-level outcome cannot be assigned to each organism as a species-specific CFU effect. Direct dental exposure, not oral supplementation. | Choi et al. 2018, [DOI 10.1016/j.bjm.2017.04.003](https://doi.org/10.1016/j.bjm.2017.04.003); [PMID 28916389 / full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC5790572/). |
| TM21 | *C. albicans* ← EDTA | `BIOFILM_PREVENTION`: concentration-dependent inhibition, XTT/microscopy, 2.5–250 mM. Preformed 24 h biofilms were minimally affected: maximum 31% XTT reduction at 250 mM. | M; 96-well in-vitro model. Mature-biofilm metabolic reduction, not proven clearance. No oral EDTA formulation or colonic delivery tested. | Ramage et al. 2007, [DOI 10.1007/s11046-007-9068-x](https://doi.org/10.1007/s11046-007-9068-x); [PMID 17909983](https://pubmed.ncbi.nlm.nih.gov/17909983/). |
| TM22 | Deoxycholate-induced *C. difficile* 630Δerm biofilms ← DNase I (proteinase K separately tested) | `DISPERSAL`: treating 48 h biofilms dispersed structure; 24 h enzyme exposure did **not** affect bacterial viability. | M; anaerobic in-vitro gut-associated chemical context. Distinct from prevention experiments. No enzyme ingestion/stability/colonic release tested. | Dubois et al. 2019, [DOI 10.1038/s41522-019-0087-4](https://doi.org/10.1038/s41522-019-0087-4); [PMID 31098293 / full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC6509328/). |
| TM23 | *E. coli* ATCC 25922 ← serrapeptase | `BIOFILM_PREVENTION`: CV IC50 14.2 ng/mL. Study additionally includes preformed-biofilm disaggregation and separate viability measurements. | M, partial; static laboratory strain model. Preserve formation, residual biomass and viability as distinct records during full extraction. No evidence here that an oral enteric enzyme product reaches intestinal biofilm at this activity. | Katsipis et al. 2025, [DOI 10.3390/microorganisms13081875](https://doi.org/10.3390/microorganisms13081875); [PMID 40871379 / full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC12388453/). |
| TM24 | *C. difficile* R20291, 1064, and 630Δerm cwp84::erm ← live *S. boulardii* **CNCM I-745** | `BIOFILM_PREVENTION`: coculture reduced biomass, viable bacterial burden and thickness; effect required direct contact and was absent with *S. cerevisiae* ATCC 9763 control. | M; lyophilized strain supplied by Biocodex, subsequently cultured; not a human gut trial or treatment of mature biofilm. Exact strain and live preparation essential. | Lacotte et al. 2022, [DOI 10.3390/microorganisms10061082](https://doi.org/10.3390/microorganisms10061082); [PMID 35744599 / full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC9227484/). |
| TM25 | *C. difficile* ← *L. reuteri* **DSM 17938 plus glycerol** | `PLANKTONIC_GROWTH` and ex-vivo community invasion prevention: reuterin-dependent inhibition; reduced invasion in antibiotic-perturbed human fecal mini-bioreactors. | M, partial; glycerol/substrate and reuterin-pathway strain context indispensable. Human-derived culture is **ex vivo**, not evidence of human oral efficacy. No mature-biofilm clearance endpoint. | Spinler et al. 2017, [DOI 10.1128/IAI.00303-17](https://doi.org/10.1128/IAI.00303-17); [PMID 28760934 / full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC5607411/). |
| TM26 | CRKP *K. pneumoniae* ZKP25 ← *L. plantarum* **LP1812** | `ANIMAL_COLONIZATION_PREVENTION`: fecal CRKP reduced below assay detection by day 5; acetate mechanism investigated. | M; antibiotic-depleted female mice, n=4/groups; 10^8 CFU by gavage for 3 days **before** oral CRKP challenge. This is prevention/anticolonization, not therapy of established human carriage. | Yan et al. 2021, [DOI 10.3389/fcimb.2021.804253](https://doi.org/10.3389/fcimb.2021.804253); [PMID 34976873](https://pubmed.ncbi.nlm.nih.gov/34976873/). |
| TM27 | *H. pylori* in asymptomatic positive adults ← nonviable *L. reuteri* **DSM 17648** | `HUMAN_LOAD_PROXY`: two-week oral preparation reduced 13C urea-breath-test value relative to placebo; coaggregation studied. | M, partial; postbiotic/inactivated cells, not live probiotic colonization. Lower UBT value is not a negative eradication test. | Holz et al. 2015, [DOI 10.1007/s12602-014-9181-3](https://doi.org/10.1007/s12602-014-9181-3); [PMID 25481036 / full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC4415890/). |
| TM28 | Gastrointestinal VRE carriage, predominantly VanB *E. faecium* ← *L. rhamnosus* **GG** yoghurt | `HUMAN_CARRIAGE_CLEARANCE`: 11/11 treatment completers versus 1/12 controls culture-negative by study completion. | M, partial; small renal-patient randomized trial, oral 100 g/day yoghurt for 4 weeks. Culture-based carriage endpoint, not biofilm measurement; attrition and recurrence/follow-up matter. | Manley et al. 2007, [DOI 10.5694/j.1326-5377.2007.tb00995.x](https://doi.org/10.5694/j.1326-5377.2007.tb00995.x); [PMID 17484706](https://pubmed.ncbi.nlm.nih.gov/17484706/). |
| TM29 | Gastrointestinal vancomycin-resistant *E. faecium* carriage ← *L. rhamnosus* **GG** capsules | `HUMAN_CARRIAGE_CLEARANCE`, **null**: no improvement in clearance or microbiome diversity. | M, partial; multicenter randomized double-blind placebo-controlled hospitalized-adult trial, 60 billion CFU/day for 4 weeks; high spontaneous clearance. Exact formulation differs from TM28. | 2022, [DOI 10.1128/spectrum.02348-21](https://doi.org/10.1128/spectrum.02348-21); [PMID 35475684 / full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC9241610/). |
| TM30 | Intestinal *Ruminococcus gnavus* signal in normal mice ← berberine | `COMMUNITY_ABUNDANCE`, downward: abundance in terminal ileum/large intestine negatively associated with increasing oral berberine exposure. | M, partial; male C57BL/6 mice, 0–300 mg/kg gavage for 2 weeks; targeted bacterial quantification. **Not MIC, direct selective killing, biofilm evidence, or human eradication.** | 2016, [DOI 10.1186/s12906-016-1367-7](https://doi.org/10.1186/s12906-016-1367-7); [PMID 27756364 / full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC5070223/). |

### 22.3 Required contradiction and preparation handling

**R. gnavus must not receive inherited antibacterial edges from E. coli, Klebsiella or Candida.** TM30 supports an animal community-abundance edge only. Direction is context-dependent: in rats with ischemia–reperfusion AKI, berberine treatment was reported to increase *R. gnavus*; add a separate upward community-abundance edge from [PMID 38291383](https://pubmed.ncbi.nlm.nih.gov/38291383/), [DOI 10.1186/s12906-023-04323-y](https://doi.org/10.1186/s12906-023-04323-y), [primary full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC10826000/). A separate metabolic-syndrome rat comparison found higher abundance in berberine than decocted *Coptis chinensis* groups ([DOI 10.1016/j.jtcms.2017.05.005](https://doi.org/10.1016/j.jtcms.2017.05.005)); the comparator must remain explicit. Neither supplies a universal beneficial/harmful interpretation of the species.

**Allicin ≠ all garlic preparations.** Human primary formulation work measured allicin bioequivalence via exhaled allyl methyl sulfide, with variable release across garlic products and meals. Enteric tablets ranged 36–104% bioequivalence with the lower-protein meal, falling to 22–57% with a high-protein meal. This verifies formulation dependence, not intact allicin concentration against gut biofilms. [Lawson & Hunsaker 2018, DOI 10.3390/nu10070812](https://doi.org/10.3390/nu10070812), [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC6073756/). Store alliin content, alliinase activity, allicin yield, processing and actual release data separately; “enteric” is a release claim, not proof of colonic activity.

**Wormwood species/preparations are not interchangeable.** Additional direct EPEC/K-12 prevention evidence exists for essential oils from *A. absinthium*, *A. campestris* and *A. herba-alba*: CV attachment reduced up to 45%/70% after 24 h at 0.62 mg/mL. [Mathlouthi et al. 2021, PMID 33588649](https://pubmed.ncbi.nlm.nih.gov/33588649/), [DOI 10.1080/08927014.2021.1886278](https://doi.org/10.1080/08927014.2021.1886278). The “up to” effects must not be assigned identically to all three oils without extracting species-specific figures. No support here for copying to *A. annua*, artemisinin, or an unspecified “wormwood” supplement.

**Mature Candida matrix caveat:** TM10 also measured Congo-red-bound matrix after treating preformed biofilm; the reduction was nonsignificant. Store an inconclusive matrix-content endpoint, not eradication.

### 22.4 Additional high-yield direct papers for next ingestion

These additional sources must be ingested at their supported extraction level. Retain incomplete fields explicitly; do not invent numerical effects.

1. *Allium ursinum* and *A. oschaninii* extracts against *K. pneumoniae* ATCC 10031 and *C. albicans* ATCC 90028 mono/polymicrobial static and dynamic biofilms: [PMID 32120894](https://pubmed.ncbi.nlm.nih.gov/32120894/), [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC7143215/). Do not call these *A. sativum* or purified allicin.
2. Methanol/ethanol *Allium sativum* extracts against six species including *E. coli* and *K. pneumoniae*, planktonic/biofilm experiments: [primary full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC4600595/). Extract organism-specific results and concentrations before representing selective effects.
3. NAC against monospecies *E. faecalis* growth and established biofilm: [PMID 22152626](https://pubmed.ncbi.nlm.nih.gov/22152626/), [DOI 10.1016/j.joen.2011.10.004](https://doi.org/10.1016/j.joen.2011.10.004); dental context.
4. *C. difficile* biofilms across five lineages, exogenous DNase, vancomycin sensitization and spore-germination measurements: [DOI 10.1038/s41598-020-78437-5](https://doi.org/10.1038/s41598-020-78437-5), [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC7865049/). Do not merge germination-based spore readout with eradication of infection.
5. *L. rhamnosus* Lcr35 inhibits EPEC, ETEC and *K. pneumoniae* adhesion to Caco-2 cells: [PMID 11316370](https://pubmed.ncbi.nlm.nih.gov/11316370/). This is epithelial competition in vitro, not durable human colonization.
6. Thymoquinone antiadhesion/antimicrobial screening including *K. pneumoniae*, antibiotic combinations: [primary full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC7989283/). Earlier broad thymoquinone prevention screen [DOI 10.1186/1472-6882-11-29](https://doi.org/10.1186/1472-6882-11-29), [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC3095572/). Preserve organisms with weak/no activity rather than projecting broad-spectrum effects.
7. Lauric acid / virgin coconut oil fractions inhibit vegetative *C. difficile*: [PMID 24328700](https://pubmed.ncbi.nlm.nih.gov/24328700/). Preparation and free-fatty-acid versus oil identity matter; no biofilm claim.

### 22.5 Database/ontology notes

Minimum evidence unit should carry: DOI/PMID/PMCID; publication status and correction/retraction check; exact preparation and compound identity; target taxonomy plus strain/pathotype; isolate body site; experiment host/site; exposure route; concentration/unit/duration; biofilm age, surface, mono/mixed culture, pH/medium; comparator; assay; result direction/value/uncertainty; viability versus biomass distinction; species-specific versus whole-community readout; human/animal/ex-vivo/in-vitro tag; and statement-level provenance. Unknown is preferable to invented defaults. Keep antibiotic/combinatorial dependencies explicit.

Useful resource specifics:

- [aBiofilm primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC5753393/), [DOI 10.1093/nar/gkx1157](https://doi.org/10.1093/nar/gkx1157), catalogs agents, target organisms, efficiency and biofilm stage. Use as discovery/index, then resolve source paper. Publication has CC BY-NC terms; a distinct dataset license and current download availability still need verification, so do not infer unrestricted commercial reuse from public accessibility.
- [BiofOmics primary paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC3386978/), [PMID 22768184](https://pubmed.ncbi.nlm.nih.gov/22768184/), is useful for assay metadata (surface, media, pH, temperature, shear), even if its live service availability needs checking. [MIABiE standards paper](https://academic.oup.com/femspd/article/70/3/250/567525) is useful for schema completeness.
- PubChem official [BioAssay documentation](https://pubchem.ncbi.nlm.nih.gov/docs/bioassays), [downloads](https://pubchem.ncbi.nlm.nih.gov/docs/downloads) and [source reuse metadata](https://pubchem.ncbi.nlm.nih.gov/docs/data-sources): record AID/SID/CID, activity outcome and assay protocol. Source-defined reuse terms mean “PubChem = blanket public domain” is unsafe.
- [ChEMBL official site](https://www.ebi.ac.uk/chembl/) states CC BY-SA 3.0 data licensing. Organism/cell assays should not be confused with biochemical protein inhibition. Preserve standard relation, units, assay type/description and source document.
- NPASS target classes include individual proteins and whole organisms; those are distinct edge types. Its 2023 paper also contains computed activity profiles/ADMET, which must remain predicted instead of measured. Use the 2026 NPASS 3.0 adapter defined in section 19.3.

### 22.6 Explicit remaining limits

Several rows are deliberately abstract-supported seeds (A), with dose/strain/assay precision still null. These are suitable for evidence discovery cards carrying their limited metadata; they are not ready for numerical cross-study ranking. Numeric efficacy rankings would additionally require comparable inocula, media, exposure, endpoint and normalization. No direct human intestinal biofilm removal evidence for oral allicin, oregano oil, wormwood, thymoquinone, serrapeptase or EDTA was established in this search. That is a search result limitation, not a proof that no such paper exists.


## 23. Butyrate producers, cross-feeding and direct precursors

BU01–BU12 extend the same candidate engine. Match exact strain and formulation before transferring any outcome; retain producer capacity, intestinal delivery and human effects separately.

### 23.1 AOR Probiotic-3: verified product identity

Official US label: https://aor.us/product/probiotic-3/ (retrieved current page). Per **two capsules**: **Enterococcus faecium T-110, 36 million CFU; Clostridium butyricum TO-A, 1.2 million CFU; Bacillus subtilis TO-A, 1.2 million CFU**. Other ingredients: lactose, potato starch, polyvinyl alcohol, polyvinylpyrrolidone; hypromellose/purified-water capsule. Contains milk (lactose). These are label composition facts, not a recommended regimen. US page says refrigerate after opening.

AOR explicitly connects its product to TOA's proprietary **Bio-Three** ingredient: https://aor.us/ingredients/bio-three-enterococcus-faecium-t-110/ . This is a manufacturer identity bridge, not a demonstration that the finished AOR capsule, its excipients, microbial proportions, or human exposure match every BIO-THREE trial tablet or veterinary formulation. Clinical BIO-THREE publications use older organism labels and frequently report microbial powder masses, which cannot simply be converted to AOR capsule CFUs.

### 23.2 Identity and archive nodes

| Node | Supported mapping | Archive / primary identity source | Boundary |
|---|---|---|---|
| C. butyricum TO-A / CBTOA | TOA strain; deposit **FERM BP-10866** | EFSA 2022: https://doi.org/10.2903/j.efsa.2022.7342 ; 2025 primary culture paper: https://pmc.ncbi.nlm.nih.gov/articles/PMC11988312/ ; genome sequences cited in that paper **CP014704, CP014705** | Distinct from CBM588 and every unspecified retail C. butyricum. Paper calls the two sequence records chromosomes; don't infer all current replicon metadata without inspecting the records. |
| E. faecium T-110, historical product label | BIO-THREE deposit **FERM BP-10867** is allocated to **E. lactis** by WGS analysis | EFSA 2022 source above; official EU identity bridge: https://joint-research-centre.ec.europa.eu/reports-and-technical-documentation/fad-2020-0058_en | Preserve label name and genome taxonomy separately. Old publications also call it Streptococcus/Enterococcus faecalis; do not describe this as a universal synonym between faecalis, faecium, and lactis. |
| T-110 published complete genome | Historical E. faecium T-110 sequence **CP006030** | Original genome paper PMID **25619602**, DOI **10.1016/j.jgg.2014.07.002**: https://pubmed.ncbi.nlm.nih.gov/25619602/ | Archival strain reference; no independently verified current AOR lot-to-genome match. |
| B. subtilis TO-A / old B. mesentericus TO-A | Deposit **FERM BP-07462**; the old name is traceable for this specific TO-A strain | EFSA 2022; Sato 2014: https://www.jstage.jst.go.jp/article/bpb/37/1/37_b13-00641/_html/-char/en | Avoid treating all historical B. mesentericus as B. subtilis. |
| AOR's three sequence links | **AB687550, AB687551, AB687552** | AOR ingredient page; AB687550 explicitly identified as 16S in Sato 2014; AB687551 NCBI title is partial 16S: https://www.ncbi.nlm.nih.gov/nuccore/AB687551 ; author preprint table labels all three 16S: https://www.preprints.org/manuscript/202309.0766 | AOR labels these “full genomic sequence data”; do **not** repeat that claim. At least AB687550/551 are directly verified gene records, not genomes. AB687552 full NCBI record was not accessible during this audit; the 16S classification is supported by an author preprint rather than independently inspected NCBI content. |
| C. butyricum MIYAIRI 588 / CBM588 | Separate strain, **FERM BP-2789** | Butirrisan study: https://pmc.ncbi.nlm.nih.gov/articles/PMC12113862/ ; official PharmExtracta portfolio: https://www.pharmextracta.com/wp-content/uploads/2023/03/PharmExtracta_Product_Portfolio.pdf | Verified Butirrisan identity; MIYA-BM also explicitly used in the 2026 study below. Neither is TO-A. |
| C. butyricum CB-a | Chicken-feces isolate; **CCTCC M 20242645** | 2026 primary paper: https://journals.asm.org/doi/10.1128/spectrum.01351-26 ; associated study BioProject **PRJNA1501883** https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1501883 | Research strain, not a verified retail product. BioProject is a study archive, not necessarily a pure-isolate complete genome. |
| WBF-011 / Pendulum Glucose Control | C. butyricum **WB-STR-0006**, A. hallii **WB-STR-0008** (old E. hallii), C. beijerinckii **WB-STR-0005**, A. muciniphila **WB-STR-0001**, B. infantis **100**, plus prebiotic | Official formulation bridge: https://pendulumlife.com/blogs/news/pendulum-therapeutics-announces-publication-of-clinical-data-for-first-ever-medical-probiotic-that-provides-the-dietary-management-of-healthy-a1c-and-blood-glucose-levels ; paper https://link.springer.com/article/10.1186/s12866-021-02415-8 ; **PRJNA722306**, GenBank **CP073277–CP073284**; metabolomics **MTBLS2713** | Archive range belongs to the study's five strains collectively; don't assign individual records without inspecting Figure 1e. Formula study, not evidence that any one component caused the outcome. |
| A. soehngenii CH-106 / CH106 | Tetracycline-sensitive mutagenized derivative of **L2-7**, distinct strain | 2025 paper https://pmc.ncbi.nlm.nih.gov/articles/PMC12087665/ ; CH106 **PRJNA1219376** (R6M8), parental L2-7 **PRJEB22345** | Do not merge with A. hallii WB-STR-0008. Developer: https://caeluslifesciences.com/ ; no verified retail product in this audit. |
| Tributyrin / CoreBiome | Defined butyrate precursor ingredient, not live bacteria | Official identity: https://compoundsolutions.com/ingredients/corebiome/ ; exact test formulation paper below | Does not require bacterial survival. Not automatically a “postbiotic” under definitions requiring inactivated microbes/components. Retail formulations need their own ingredient and release-system match. |

### 23.3 Twelve exact evidence seeds

#### BU01 — TO-A lactate/acetate conversion, 2025

- **Primary:** Honda et al.; PMID **40243571**; DOI **10.3390/ijms26072951**. https://pmc.ncbi.nlm.nih.gov/articles/PMC11988312/
- **Intervention / model:** Exact C. butyricum TO-A in anaerobic peptone–yeast cultures, with defined organic-acid substrates; three replicates.
- **Target / actual effect:** Directly measured butyrate production. Net production after 16 hours was **2.59 ± 0.10 mM** in base medium, **4.94 ± 0.83 mM** with acetate plus D/L-lactate, and **7.54 ± 0.67 mM** with acetate plus D-lactate. Supports a cross-feeding mechanism using lactate and acetate; relevant genes and racemization were investigated.
- **Gut delivery:** None measured; cultured viable counts are not recovery from a gastrointestinal tract.
- **Limits:** Manufacturer-affiliated authors; simplified culture. Does not establish human colonic delivery, a fecal increase, or efficacy of finished AOR capsules. Avoid interpreting the abstract's “only” as inability to make butyrate otherwise—the control medium produced butyrate too.

#### BU02 — Veterinary BIO-THREE increases measured fecal butyrate, 2020

- **Primary:** Inatomi, Makita, Otomaru; DOI **10.12935/jvma.73.374**. https://www.jstage.jst.go.jp/article/jvma/73/7/73_374/_article/-char/en ; full primary PDF https://www.jstage.jst.go.jp/article/jvma/73/7/73_374/_pdf/-char/en
- **Intervention / model:** Forty healthy eight-day-old Japanese Black calves, randomized 20/20; veterinary BIO-THREE from TOA identified in Japanese methods; 14-day experiment.
- **Target / actual effect:** Fecal GC assay. Butyrate **20.7 ± 2.7 versus 16.4 ± 3.3 µmol/g** on day 7 and **21.6 ± 3.7 versus 17.3 ± 2.9** on day 14; both p<0.05. Fecal IgA also increased; lactate did not differ significantly.
- **Gut delivery:** Fecal metabolite result, no strain-specific viable recovery assay.
- **Limits:** Veterinary mixture, suckling calves, single farm; no human translation of magnitude. Primary paper names species/product but not strain suffixes, so connect through the documented TOA formulation and retain formulation-version uncertainty. Does not isolate TO-A's contribution.

#### BU03 — BIO-THREE ulcerative-colitis remission RCT, 2015

- **Primary:** Yoshimatsu et al.; PMID **26019464**; DOI **10.3748/wjg.v21.i19.5985**. https://www.wjgnet.com/1007-9327/full/v21/i19/5985
- **Intervention / model:** Exact BIO-THREE tablets listing T-110, C. butyricum TO-A and B. mesentericus TO-A; 60 randomized patients with inactive UC, 46 completed; adjunct to ongoing medicines.
- **Target / actual effect:** At 12 months remission was **69.5% versus 56.6%, p=0.248**. The early three-month relapse contrast favored probiotic; later comparisons were nonsignificant. Fecal SCFA ratios and bacterial profiles were examined.
- **Butyrate interpretation:** This is **not** evidence that increasing stool butyrate mediates benefit. Higher fecal butyrate/acetate ratios preceded relapse; authors discuss impaired absorption/utilization. Their earlier active-UC findings associated improvement with *lower* fecal butyrate.
- **Gut delivery:** Community DNA profiling, not culture-confirmed survival of each administered strain.
- **Limits:** Small completed-case sample, attrition, exploratory tiny microbial subgroups; main long-term result nonsignificant. Finished AOR-capsule equivalence untested.

#### BU04 — BIO-THREE dialysis RCT, 2025

- **Primary:** Ogawa et al.; DOI **10.1186/s41100-025-00657-0**. https://link.springer.com/article/10.1186/s41100-025-00657-0
- **Intervention / model:** Named BIO-THREE three-strain product; maintenance-dialysis patients; 73 analyzed (37 probiotic, 36 placebo), six months, blinded randomized design.
- **Target / actual effect:** **No significant change in primary hs-CRP/albumin endpoint.** Serum phosphorus declined within the probiotic group; Lachnoclostridium abundance increased and correlated with phosphorus changes. Laxative-related findings were exploratory.
- **Butyrate:** No direct butyrate assay supporting an increase. Calling a taxon “butyric acid bacteria” does not substitute for metabolite measurement.
- **Gut delivery:** 16S community sequencing; no administered-strain viable recovery.
- **Limits:** Modest size, baseline imbalances and incomplete adherence; within-group phosphorus change is weaker than a clear between-group treatment effect. Relevant exact-mixture human evidence, not proof of systemic anti-inflammatory benefit.

#### BU05 — BIO-THREE / MIYA-BM COVID observational update, 2026

- **Primary:** Morikawa et al.; PMID **42127027**; DOI **10.1080/21505594.2026.2673650**. https://pubmed.ncbi.nlm.nih.gov/42127027/ ; https://www.tandfonline.com/doi/full/10.1080/21505594.2026.2673650 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC13182968/
- **Intervention / model:** Retrospective hospital records from 2020–2021, reported May 2026; 27 C. butyricum-product recipients versus 256 nonrecipients. BIO-THREE subgroup **n=10**; MIYA-BM is a separate preparation.
- **Target / actual effect:** Mortality, ventilation duration, secondary bacterial pneumonia. BIO-THREE subgroup showed favorable numerical differences but **no significant differences** for these outcomes (e.g. mortality 0 versus 15.6%).
- **Butyrate / delivery:** Authors explicitly did not measure fecal C. butyricum, SCFAs, or microbiome profiles.
- **Limits:** Very small subgroups, treatment selection/confounding, old clinical period. Keep as a current exact-product observational seed; no antiviral, butyrate-production, or superiority claim is established.

#### BU06 — Exact B. TO-A human viable recovery, 2014

- **Primary:** Sato, Seo, Benno; DOI **10.1248/bpb.b13-00641**. https://www.jstage.jst.go.jp/article/bpb/37/1/37_b13-00641/_html/-char/en
- **Intervention / model:** Standalone **B. mesentericus TO-A** tablets, 24 healthy volunteers; seven-day administration and seven-day withdrawal. The study explicitly excluded other BIO-THREE products during the experiment.
- **Target / actual effect:** Fecal strain detection by culture plus strain-specific qPCR. All participants positive during administration; averages approximately **10^4.8 CFU/g by culture** and **10^5.8 cells/g by qPCR**. Half remained qPCR-positive three days after withdrawal; only one remained positive after one week.
- **Gut delivery:** **Direct human viable recovery for this TO-A Bacillus component**, stronger than DNA-only detection. Supports transient passage, not persistent engraftment.
- **Butyrate:** Not measured and Bacillus recovery does not establish C. butyricum TO-A delivery.
- **Limits:** Different standalone product and study exposure; not a test of current AOR capsules. Manufacturer involvement.

#### BU07 — CBM588 rat germination / intestinal fate, 1997

- **Primary:** Sato and Tanaka; PMID **9343816**; DOI **10.1111/j.1348-0421.1997.tb01909.x**. https://onlinelibrary.wiley.com/doi/full/10.1111/j.1348-0421.1997.tb01909.x
- **Intervention / model:** Oral **C. butyricum MIYAIRI 588 (CBM588)** spores in rats, selective culture, ethanol-treated spore assay and specific antibody imaging.
- **Target / actual effect:** Evidence of germination in upper small intestine, vegetative cells distally after two hours, growth in cecum/colon after five hours; no intestinal detection after three days.
- **Gut delivery:** **Direct animal viability and localization**, strain-specific.
- **Butyrate:** No quantitative butyrate endpoint in the accessed abstract.
- **Limits:** Animal physiology and single administration; transient fate. **Do not attach this as direct TO-A survival evidence.** The 2025 TO-A EFSA dossier mentions it but explicitly identifies it as a different strain.

#### BU08 — CB-a barrier / SCFA study, 2026

- **Primary:** Liu et al.; PMID **42549889**; DOI **10.1128/spectrum.01351-26**; published August 4, 2026. https://journals.asm.org/doi/10.1128/spectrum.01351-26
- **Intervention / model:** Exact **C. butyricum CB-a**, CCTCC M 20242645; 40 male mice in four groups, DSS injury model.
- **Target / actual effect:** Measured fecal SCFA restoration, particularly butyrate (p<0.05), alongside reduced inflammatory markers and restoration of ZO-1/occludin, altered GPR41/43/109A and HDAC1/2 measures. This supports a barrier-and-metabolite mechanism to investigate.
- **Gut delivery:** Oral administration plus community/metabolite observations; not a demonstrated human viable-recovery result.
- **Limits:** Mouse study, strain isolated from chicken feces, no verified retail formulation. Correlated pathway measurements do not by themselves prove butyrate is the sole necessary mediator. No extrapolation to TO-A, CBM588, or unspecified C. butyricum supplements.

#### BU09 — WBF-011 / Pendulum actual human plasma butyrate, 2022

- **Primary:** McMurdie et al.; DOI **10.1186/s12866-021-02415-8**. https://link.springer.com/article/10.1186/s12866-021-02415-8
- **Intervention / model:** WBF-011 five-strain plus prebiotic formulation in a 12-week randomized placebo-controlled T2D trial; 76 randomized, 51 with paired plasma suitable for metabolomics.
- **Target / actual effect:** Fasting **plasma butyrate median increase 0.15 µM (~27%)**, nominal one-sided between-group p<0.05. This is a human measured-metabolite result, not merely an inferred producer abundance. The original fecal increase was nonsignificant.
- **Gut delivery:** Systemic metabolite exposure supports metabolic activity but does not quantify colon production, anatomical delivery, or viable recovery of every strain.
- **Limits:** Exploratory analysis and mixture with prebiotic; attribution to individual strains is impossible. Drug-associated effect heterogeneity and manufacturer authorship matter. Do not transfer findings to a single-strain Pendulum product. Study archives: PRJNA722306, CP073277–CP073284, MTBLS2713.

#### BU10 — A. soehngenii L2-7 human trial with SCFA null, 2024

- **Primary:** Attaye et al. https://pmc.ncbi.nlm.nih.gov/articles/PMC11321313/
- **Intervention / model:** Exact **A. soehngenii L2-7**; 25 men with T2D on stable metformin; 14-day blinded randomized placebo-controlled study.
- **Target / actual effect:** Improved glycemic-variability measures; **no notable changes in measured SCFAs**. A valuable producer candidate with an explicitly retained human metabolite null result.
- **Gut delivery:** Oral intervention and microbiome measurements; do not convert fecal DNA detection into culture-proven viability.
- **Limits:** Small short trial in men receiving metformin. Clinical outcomes cannot be assumed to be caused by increased butyrate. No verified off-the-shelf retail product. Keep distinct from CH106 and A. hallii WB-STR-0008.

#### BU11 — A. soehngenii CH106 next-generation strain, 2025

- **Primary:** DOI **10.1080/19490976.2025.2504115**, PMID **40371708**. https://pmc.ncbi.nlm.nih.gov/articles/PMC12087665/
- **Intervention / model:** Encapsulated freeze-dried **CH-106**, 98 insulin-resistant/prediabetic adults, three-month blinded randomized placebo-controlled multicenter trial.
- **Target / actual effect:** Glycemic variability and diastolic pressure improved; HbA1c difference emerged after 16 weeks including follow-up. Exact strain produces butyrate in culture from sugars or lactate plus acetate; no direct demonstrated human butyrate increase in this trial.
- **Gut delivery:** Fecal strain-specific qPCR increased approximately 100-fold. Active-cell formulation checks support product viability before ingestion; fecal qPCR remains **DNA recovery**, not culture-confirmed survival.
- **Limits:** Geographic response heterogeneity, modest clinical effect sizes, multiple mechanisms plausible. Manufacturer/developer involvement. CH106 has 37 genome mutations relative to L2-7, including tetO inactivation, so they need separate nodes. Archives PRJNA1219376 and parental PRJEB22345.

#### BU12 — CoreBiome tributyrin simulated delivery and butyrate, 2025

- **Primary:** Duysburgh et al.; DOI **10.3389/fnut.2025.1712993**. https://www.frontiersin.org/journals/nutrition/articles/10.3389/fnut.2025.1712993/full
- **Intervention / model:** Tested CoreBiome capsule/softgel formulations; simulated upper digestion, SHIME reactors from three healthy male fecal donors, Caco-2/THP1 assays.
- **Target / actual effect:** Reactor butyrate increased. Barrier and immune readouts changed in cell assays exposed to reactor supernatants.
- **Delivery:** At simulated upper-GI completion, **59.1% of capsule tributyrin / 51.3% of softgel tributyrin remained**; the rest was hydrolyzed in the simulated small intestine. Investigators then added liquid tributyrin directly to colon reactors using the estimated surviving fraction.
- **Limits:** **In vitro delivery estimates**, not human colonic bioavailability. Three donors, no placebo-controlled human efficacy from this paper, manufacturer funding. Useful direct precursor alternative to live organisms; release-system equivalence must be verified for any finished consumer product.

### 23.4 Delivery evidence must remain a separate field

| Question | Supported answer from these sources |
|---|---|
| Does exact C. butyricum TO-A make butyrate? | Yes, directly in culture (BU01); BIO-THREE also increases measured fecal butyrate in calves (BU02). |
| Do current AOR capsules reliably increase human butyrate? | Not demonstrated by the retrieved exact-product human sources. Keep the option with the evidence level shown. |
| Is there human live recovery for a component? | Yes, standalone Bacillus TO-A (BU06); this does not prove all three survive the finished product. |
| Is there direct animal TO-A recovery? | EFSA's 2022 chicken dossier recovered culture colonies matching C. butyricum FERM BP-10866 by RAPD in **11 of 83 sampled specimens**, but counts were low and the panel could not draw compatibility conclusions. Bacillus strain matching succeeded; Enterococcus matching did not. This is limited animal evidence, not no evidence. |
| Is the classic rat germination paper TO-A? | No: CBM588 (BU07). |
| Does fecal PCR prove a live strain reached the colon? | No; record DNA recovery separately from culture recovery and metabolic activity. |
| Does more fecal butyrate always mean better availability? | No; stool reflects production, absorption, utilization, transit and dilution. BU03 shows why a high stool ratio can coexist with worse disease. |
| Is a measured plasma increase available for another product? | Yes, WBF-011/Pendulum mixture (BU09), with exploratory-analysis limitations. |
| Can direct butyrate precursors bypass the survival question? | Yes mechanistically; BU12 adds formulation-specific simulated-release evidence, without claiming human colon-delivery percentages. |

EFSA source for the limited chicken-recovery row and taxonomy: https://efsa.onlinelibrary.wiley.com/doi/10.2903/j.efsa.2022.7342 . Treat the official dossier as an assessment of submitted experiments, not an independent human efficacy study.

### 23.5 Additional follow-up leads (not merged into the 12 seeds)

- **Butirrisan/CBM588 IBS-D, 2025:** https://pmc.ncbi.nlm.nih.gov/articles/PMC12113862/ . Open-label single-arm, 205 adults plus dietary change; symptom improvement and 19-person microbiota subset; no direct butyrate assay. Product identity verified by official PharmExtracta portfolio. Relevant option, stronger identity than an unspecified-strain marketplace supplement.
- **CBM588 metabolic engineering, 2025:** PMID39892707, DOI10.1016/j.anaerobe.2025.102940, https://pubmed.ncbi.nlm.nih.gov/39892707/ . Culture manipulation of Ptb/Buk/Crt pathways; enhanced engineered derivatives must not be presented as ordinary retail CBM588.
- **BIO-THREE oxaliplatin injury, 2021 online/2022 issue:** PMID33956306, DOI10.1007/s12602-021-09795-3, https://pubmed.ncbi.nlm.nih.gov/33956306/ . Primary human/animal/translational lead; not fully extracted in this audit.
- **TO-A/BIO-THREE renal-injury rat experiment, 2023:** DOI10.1371/journal.pone.0281745, https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0281745 . Exact mixture identified; retain as mechanistic lead rather than a butyrate-assay claim.
- **Safety belongs at strain/population level:** The 2024 primary bacteremia report genetically links some clinical C. butyricum isolates to probiotic use: https://wwwnc.cdc.gov/eid/article/30/4/23-1633_article . Do not turn this into a blanket statement that all C. butyricum strains/products are pathogenic, or into a blanket assumption of safety from spore status.

### 23.6 Suggested inclusion labels

AOR/BIO-THREE: **“Exact strain producer; human outcome studies; human butyrate increase unconfirmed.”** CBM588: **“Distinct spore-forming producer with direct animal gut-fate evidence.”** WBF-011: **“Exact mixture, measured human plasma-butyrate increase.”** L2-7: **“Next-generation producer, human SCFA null retained.”** CH106: **“Related distinct strain, human DNA recovery and clinical signals.”** CoreBiome: **“Direct precursor, simulated delivery and reactor butyrate.”** CB-a: **“Research strain, animal barrier/metabolite mechanism.”**

Avoid labeling all of these “clinically proven butyrate boosters.” Do include all as legitimate research options with explicit evidence type and identity confidence.

