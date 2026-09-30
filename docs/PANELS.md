# Pathways

Reference for the 25 shipped gene panels, how to add another, and what is
deliberately left out.

---

## The panels

Grouped the way the report groups them. "Headline from" is the
`aggregate_from` field: every gene is always reported, but the panel-level
figure is computed only from the genes whose reference sets are specific
enough to carry it.

**Metabolite production capacity**

| Panel | Metabolite | Genes | Headline from | Notes |
|---|---|---|---|---|
| `butyrate` | Butyrate | `but`, `buk` | both | Positive control — abundant in any sample |
| `propionate` | Propionate | `pct`, `pduP`, `lcdA` | `pct`, `pduP` | Succinate route excluded |
| `bcaa` | Branched-chain amino acids | `ilvC`, `ilvD`, `leuA` | `ilvC`, `ilvD` | Insulin-resistance link (Pedersen 2016) |
| `cutc` | Trimethylamine / TMAO | `cutC`, `cntA` | both | Two independent substrate routes |
| `urda` | Imidazole propionate | `urdA` | `urdA` | Has an active-site residue check |
| `bai` | Secondary bile acids | `baiCD`, `baiH`, `baiB`, `baiF`, `baiE`, `baiA` | `baiCD`, `baiH` | Rare organisms; smallest reference sets |
| `bsh` | Unconjugated bile acids | `bsh` | `bsh` | High prevalence, low discriminating power |
| `indole` | Indole / indoxyl sulfate | `tnaA` | `tnaA` | Common in Enterobacteriaceae and Bacteroidales |
| `ipa` | Indole-3-propionate | `fldBC`, `fldH` | `fldBC` | Narrow group of Clostridia |
| `pcresol` | p-Cresol / p-cresyl sulfate | `hpdB` | `hpdB` | Large subunit only |

**Neuroactive compounds**

| Panel | Metabolite | Genes | Headline from | Notes |
|---|---|---|---|---|
| `gaba` | GABA | `gadB`, `gabT` | `gadB` | Production and degradation both shown; see below |
| `dopamine` | Dopamine (from L-DOPA), tyramine | `tyrDC`, `dadh` | `tyrDC` | `dadh` has two reference sequences — read as presence only |
| `histamine` | Histamine | `hdcA`, `hdcP` | `hdcA` | Pyruvoyl HDC; the PLP-type is a decoy |
| `tryptamine` | Tryptamine | `tdc` | `tdc` | Two known gut producers |

Serotonin has no bacterial synthesis gene to count: the gut's serotonin is
host-made, and bacteria influence it indirectly (SCFAs, secondary bile acids,
tryptamine). The panels above cover those levers; a "serotonin capacity"
number would be invented.

**Vitamin synthesis**

| Panel | Vitamin | Genes | Headline from |
|---|---|---|---|
| `b12` | B12 (cobalamin) | `cbiA`, `cobS`, `btuB` | `cbiA`, `cobS` |
| `k2` | K2 (menaquinone) | `menA`, `menD`, `mqnE` | all |
| `folate` | B9 (folate) | `folP`, `folK`, `folC` | `folP`, `folK` |
| `riboflavin` | B2 | `ribB`, `ribE`, `ribH` | all |
| `biotin` | B7 | `bioB`, `bioF`, `bioD` | all |

**Gas production**

| Panel | Gas | Genes | Headline from | Notes |
|---|---|---|---|---|
| `h2s` | Hydrogen sulfide | `dsrA`, `asrA`, `cgl` | `dsrA`, `asrA` | Cysteine route reported, not in the headline |
| `methane` | Methane | `mcrA`, `mcrB`, `mtaB` | `mcrA`, `mcrB` | Archaeal; taxonomy filter is Archaea, not Bacteria |

**Detoxification and breakdown**

| Panel | Substrate | Genes | Headline from | Notes |
|---|---|---|---|---|
| `bglucuronidase` | Glucuronides (oestrogens, drug metabolites) | `gus` | `gus` | Loop-1 GUS cannot be separated at read level; whole family counted |
| `oxalate` | Oxalate | `oxc`, `frc` | both | Degradation, not production |
| `urease` | Urea → ammonia | `ureC`, `ureB` | `ureC` | Alpha subunit carries the headline |

**Toxins and virulence**

| Panel | Factors | Genes | Headline from | Notes |
|---|---|---|---|---|
| `crc_virulence` | FadA, *B. fragilis* toxin, colibactin | `fadA`, `bft`, `clbB` | all | Genus-restricted UniProt filters; tiny reference sets by nature |

Inspect any panel in full, including every UniProt query:

```bash
openbiota panels                  # summary table
openbiota panels --show urda      # one panel, complete
openbiota panels --json           # machine-readable
```

---

## Deliberately excluded

**`ldh` (lactate), `speE` (spermidine), `murI` (D-glutamate).** These genes are
near-universal across bacteria. Essentially every genome scores positive, so
the number measures how much bacterial DNA is in the sample — which `rpoB`
already reports. A test asserts they stay absent.

**`gadB` was on that list and is now screened, with a condition.** Glutamate
decarboxylase is widespread but not universal — it is concentrated in
Bacteroides, Bifidobacterium, Lactobacillus and Enterobacteriaceae and absent
from most Clostridia — and its spread across the reference cohort is what
decides whether the number carries information. The panel reports the
degradation gene `gabT` beside it because the net GABA available depends on the
balance. If a cohort rebuild ever shows the panel's coefficient of variation
collapsing towards the normaliser's, it goes back on the excluded list.

**Urolithin A cannot be screened at all.** The bacteria that convert
ellagitannins to urolithin A have not been identified, so there is no marker
gene. The `polyphenol_metabolisers` taxon group reports the known converter
species (*Gordonibacter*, *Ellagibacter*) instead, as carriage rather than
capacity.

**The succinate route to propionate.** Its marker genes sit in the core
methylmalonyl-CoA pathway, so nearly every genome would score positive. The
`propionate` panel therefore covers the acrylate and propanediol routes only,
and understates total propionate capacity.

**Serotonin, melatonin, acetylcholine.** No bacterial biosynthesis gene that is
both specific and present in gut genomes. Reporting a number would be
reporting a guess.

---

## Adding a pathway

A pathway is a YAML file in `panels/`. No code changes.

```yaml
name: mypanel
description: What this pathway does.
metabolite: Short label            # fits ~33 characters in the summary table
category: metabolite               # metabolite | neuroactive | vitamin | gas | detox | virulence

targets:
  - id: MYGENE
    label: myGene (full protein name)
    gene: myGene
    source: uniprot
    query: '(protein_name:"my protein name") AND (taxonomy_id:2)'
    min_identity: 60.0

decoys:
  - id: HOMOLOG
    label: the close homolog that would otherwise produce false positives
    source: uniprot
    query: '(protein_name:"homolog name") AND (taxonomy_id:2)'
    min_identity: 50.0

interpretation:
  what_it_is: Plain-language description for the consumer report.
  made_from: Dietary precursor
  higher_means: adverse            # adverse | favourable | context-dependent | unclear
  evidence_strength: limited       # strong | moderate | limited | preliminary
  summary: One or two sentences on what the research shows.
  evidence_detail: The specific findings, with study context.
  citation: Author et al., Journal (Year), doi:...

citation: >-
  Author et al., Journal (Year), doi:...
```

Then:

```bash
openbiota panels --show mypanel     # validate and inspect
openbiota build-db                  # fetch and check the reference sets
openbiota run --subsample 100000    # smoke test
openbiota validate                  # confirm it does not degrade accuracy
```

### Optional keys

| Key | Default | Purpose |
|---|---|---|
| `pathway` | — | Reaction sequence, shown in the summary |
| `organisms` | — | Literature-reported producers |
| `aggregate` | `sum` | `sum`, `median`, `max`, `mean` |
| `aggregate_from` | all targets | Narrow the headline to your specific genes |
| `aggregate_reason` | — | Why that aggregate — appears in the report |
| `min_fragments_for_stability` | 20 | Below this, flag as statistically unstable |
| `expected_copies_per_genome` | 1 | Typical copies of one headline gene per carrier genome; raises the `implausible-*` ceiling for paralog-rich families such as GUS |
| `min_alignment_aa` | 25 | Minimum aligned residues per read |
| `length_range` | — | Explicit `[min, max]` reference length window |
| `length_tolerance` | 0.25 | Median-relative window when `length_range` is absent |
| `max_sequences` | 4000 | Cap on references fetched per entry |
| `residue_check` | — | Active-site residue filter (see below) |
| `caveats` | — | Panel-specific notes, printed in the summary |
| `extension` | `false` | Marks a panel as beyond the original scope |

### Residue check

```yaml
residue_check:
  enabled: true
  anchor_accession: Q8CVD0       # a well-characterised reference protein
  canonical_position: 373        # position in the anchor
  accepted_residues: [Y, M]
  applies_to: [URDA]
  description: What this residue does.
```

At build time every reference is aligned to the anchor, and its own equivalent
position and residue are recorded. The build aborts if the anchor fails to map
its own canonical position onto itself.

---

## Four rules for writing queries

Learned from the queries that went wrong.

**1. Query on the recommended protein name, not the gene name.** UniProt gene
symbols are heavily reused across unrelated proteins. `gene:cutC` returns ~6,300
entries that are mostly the *copper homeostasis* protein CutC; `gene:cntA`
returns a metal-staphylopine transporter. The recommended-name assignment is a
curation layer built on family membership and signature residues. Gene symbols
are not.

**2. Never put a comma inside parentheses inside a quoted phrase.** It silently
breaks UniProt's boolean parser and your `taxonomy_id` filter stops applying.
You get a smaller, plausible, wrong number and no error at all. A test enforces
this across every shipped query.

**3. Check the organisms your query returns, not just the count.** The obvious
`baiCD` query returned 361 hits that were mostly *plant* progesterone
5β-reductases. `openbiota build-db` prints the kept count, mean length and length
window for every entry; the per-panel hit files show which organisms actually
matched. A reference set full of organisms that cannot perform the chemistry is
a reference set that will mislead.

**4. Declare the close homolog as a decoy.** Best-hit competition is the primary
specificity mechanism. An identity floor alone is not enough.

A zero-sequence reference fetch is a hard error, never a count of zero.

---

## The built-in normaliser

`rpoB` is not a panel. It lives in `panels/_normalizer.yaml` and is included in
every run.

```yaml
normalizer:
  id: RPOB
  query: '(gene:rpoB) AND (taxonomy_id:2) AND (reviewed:true)'
  min_identity: 50.0
  length_range: [1000, 1500]
```

Reviewed-only, because a fragmentary or misannotated reference would distort the
mean reference length in the denominator, and reviewed `rpoB` already spans
every major bacterial phylum. The explicit length window also rules out `rpoC`
(subunit beta-prime) leaking in through name similarity.
