# Expanded gut bacterial and archaeal recognition (spec 0.8.4)

This document says what the detection expansion is made of, where every
piece came from, how it is installed, and what each lane may and may not
claim. The numbers a run actually used are in `results.json` under
`detection.reference_releases`, in `refs/expanded/locks/*.lock.json`, and in
`refs/expanded/tools.lock.json`.

## Install

```
make tools-expanded    # every pinned tool, then scripts/lock_tools.py
make refs-expanded     # every reference, resumable and checksummed, then the Jan26 index
make expansion         # crosswalks + supplement reconciliation + rescue panel manifest
make refs-globdb-genomes  # the 272 GB GlobDB genome archive, fetched, verified and unpacked (one-time)
make expansion-cohorts    # a reference population per expanded lane, for percentiles in that lane's namespace
```

`scripts/fetch_expanded_refs.py --status` shows download progress. A file
whose lock verifies is never downloaded again. Files over 1 GB are fetched
with `aria2c` (eight parallel range requests, resumable; the 292 GB GlobDB
archive finishes in minutes rather than hours) and every file's SHA-256
and MD5 are computed in one pass when it lands.

`scripts/refs_manifest.py` writes `refs/MANIFEST.json` and `refs/MANIFEST.md`:
one row per directory under `refs/` with its size, file count, category
(locked download, derived index, local cohort, on-demand cache), the lock
that covers it and the command that rebuilds it. `--check` exits non-zero
if any lock-covered file is missing or the wrong size. This is the packing
list for archiving the repository or installing it on another machine.

## Tools (pinned; `refs/expanded/tools.lock.json`)

| Tool | Version | Installed by | Used for |
|---|---|---|---|
| MetaPhlAn | 3.1.0 | `.venv-mpa3` | scoring lane (every percentile) — unchanged |
| MetaPhlAn | 4.1.1 | `.venv-mpa4` | Jun23 baseline lane — unchanged |
| MetaPhlAn | 4.2.6 (tag c55b299; reports 4.2.5) | `.venv-mpa42` | Jan26 marker lane, sample2markers, StrainPhlAn |
| sylph | 1.0.0 | `vendor/sylph` | GTDB R232 lane (unchanged) and GlobDB r232 lane |
| mOTUs | 4.1.0 | `.venv-motus` (python 3.12) | universal-marker lane |
| Kraken2 / Bracken | 2.1.x / 3.x | brew | gut rescue classifier |
| SingleM | 0.21.4 | `.venv-singlem` (+ orfm, smafa, mfqe) | unrepresented-lineage rescue |
| skani, fastANI | 0.3.x, 1.3x | brew | genome ANI + alignment fraction |
| bowtie2, samtools | 2.5.x, 1.x | brew | competitive confirmation |
| mafft, trimal, RAxML, blastn | brew / vendor | via `vendor/bin` | PhyloPhlAn inside StrainPhlAn |
| SPAdes | 4.3.0 (Darwin binary) | `vendor/` | targeted assembly (metaSPAdes) |

## References (locked; `refs/expanded/locks`)

| ID | Source | Release | Role |
|---|---|---|---|
| R1 | MetaPhlAn CHOCOPhlAnSGB | `mpa_vJan26_CHOCOPhlAnSGB_202605` (72,000 SGBs) | marker lane |
| R2 | GTDB | R232 (199,923 species clusters) | backbone; every canonical species name |
| R3 | GlobDB | r232 (346,233 representatives, 26 sources) | broad sylph discovery; source dictionaries |
| R4 | mOTUs DB | 4.1 (124,295 mOTUs in the artifact; GTDB R226 names) | universal-marker lane |
| R5 | HRGM2 | 4,824 species representatives | supplement; benchmark truth genomes |
| R6 | UHGG | v2.0.2 (4,744 species; publisher Kraken2 DB) | gut rescue classifier |
| R7 | ELGG | 2,172 representatives | supplement (early-life MAGs) |
| R8 | HumGut2 | 31,225 genomes | supplement (metadata; genomes on demand) |
| R9 | HROM | 5,113 oral species | supplement (metadata; genomes on demand) |
| R10 | NCBI | accession-versioned | genomes for confirmation and strain references |
| — | SingleM metapackage | GlobDB_r232.metapackage_v4 | lineage rescue |
| — | Upstream SGB bridges | Jun23→R207, Jan26→R226 (Jan25→R220 kept only for audit) | crosswalks |

## Crosswalks (`refs/crosswalk`)

A GTDB species name is not stable across releases, so the upstream bridge
is only the first step. `openbiota.expansion.crosswalk` walks each
bridge species to R232 through its member genomes and types the result:
`exact`, `synonym` (renamed), `split` (rendered as a complex, never as
several species), `merge` (several SGBs, one organism), `parent`,
`unresolved`. Jan26→R232: 50,648 exact, 3,320 renamed, 37 splits, 74
merges. Jun23→R232: 6,137 renamed — the scale of what the misapplied R220
file had been getting wrong. mOTUs (R226 names) use the same walk.

Two more bridges keep one organism from appearing twice
(`openbiota.expansion.names`):

* **MetaPhlAn 3 names → GTDB R232 species.** The scoring lane carries 2019
  NCBI names (`Blautia_sp_CAG_257`, `Firmicutes_bacterium_CAG_424`). Each
  is given the R232 species that the genomes of that NCBI organism name
  were placed in (majority over the release's genomes, from the R232
  metadata; cached to `refs/crosswalk/ncbi_names_to_gtdb_r232.tsv.gz`), so
  the record merges with the newer catalogues' record of the same organism
  and keeps its scoring percentile through the rename. A name whose genomes
  scatter across GTDB species (`Faecalibacterium prausnitzii`) maps to
  nothing and stands under its own label.
* **mOTUs clusters → representative genomes.** The profiler's label is an
  80 % majority vote over a cluster's members and reads `Unknown <genus>`
  when they disagree. The cluster's representative genome
  (`mOTUsv4.1.gtdb.taxonomy.rep.tsv.gz`) has a GTDB species, and GlobDB
  may hold that very genome (`MOTU4_dictionary.tsv`): an unknown-labelled
  cluster is named by its representative's species when that is a species,
  otherwise by its GlobDB genome id, so the mOTUs lane and the sylph lane
  produce one record for one unit. Records that name the same reference
  genome, under any label, are folded into one.

## Reconciliation (`refs/expanded/reconciliation`)

Every supplement representative is placed in one category with the
method that put it there: `already_represented` (accession, assembly walk,
name match or ANI), `new_cluster` (a GlobDB non-GTDB representative from
that source), `additional_strain_reference`, `better_reference`,
`unresolved_mapping`, `excluded_quality`. Counts are per source and never
summed. The ANI rule for a provisional cluster is 95% ANI and 65% aligned
fraction of the shorter genome; boundary cases are flagged, not forced.

UHGG v2.0.2 representatives that no accession or name could place (470
MGnify MAGs GTDB never took in) are placed by sequence
(`scripts/reconcile_uhgg_ani.py`): each genome is queried against the
GlobDB r232 sylph database and the best hit's naive ANI and k-mer
containment decide - same species at ≥ 95.5% ANI with a quarter of the
genome contained, boundary hits (94–95.5%) recorded and left unresolved.
Without this, on the simulated communities those representatives
surfaced as "new clusters" that were truth species under another
catalogue's identifier, and in the reports as unnamed clusters beside the
named species they are.

## How a run uses the lanes

Every lane runs on the same host-filtered reads and records its own
observations in its own units (`detection.lanes`). The inventory merges
identities across lanes (aliases, GTDB names, crosswalk merges), takes the
composition share from one primary lane (Jan26), keeps every other lane's
reading as a secondary value in its own unit, and never sums across lanes.

Two Kraken lanes run the same tool at the same operating point and count
as one method: the publisher's UHGG v2.0.2 database, and the gut rescue
panel built locally by `scripts/build_rescue_panel_kraken.py` from the
reconciliation targets (R5–R10 clusters absent from UHGG), every GlobDB
gut-source representative in those targets' genera as near neighbours,
every reference genome any sample's confirmation has fetched as a detected
alternative, and PhiX/food decoys (first 200 Mb of each background genome).
Decoy hits stay in the lane record and are never organisms. The panel's
taxonomy is custom and reversible (`panel_manifest.tsv`: taxid ↔ source id
↔ lineage ↔ role).

### One share of the whole per organism (`inventory.unify_shares`)

The primary lane measures the composition but names populations by its
own catalogue, so an organism it has one entry for may be two populations
the whole-genome catalogues tell apart (*Lachnospira eligens* and *L.
eligens_A*), and a cluster it has no entry for it cannot place at all. In
the six project samples, the whole-genome lanes' readings for organisms the
primary lane did not name summed to 50–65% of their own scale while the
primary lane's unclassified band was 9–15%: almost all of that mass is
populations the primary lane had already counted under a relative's name.
After competitive confirmation the inventory therefore resolves every
organism onto the primary scale:

* **split** — a marker species is divided among the populations the marker
  catalogue cannot tell apart from it: its GTDB siblings (*X*, *X_A*, *X_B*)
  that the competition distinguished, in proportion to the reads that
  mapped uniquely to each genome per megabase. The species total is
  conserved and the lane's own reading is kept as `marker_percent`; the
  marker lane's readings for species it does distinguish are never
  redistributed between them.
* **estimated** — an organism whose genus the primary lane did not place at
  all takes its whole-genome reading rescaled to the primary scale (median
  ratio of the two lanes over organisms both measured, 0.75–0.8 here);
  those shares are drawn from the unclassified band (`unplaced_percent` is
  what remains).
* **member** — a population the whole-genome lanes saw inside a genus the
  primary lane measured that is not a sibling of any marker species there
  (an unnamed cluster, a species the marker catalogue lacks), or that the
  competition could not apportion: its reads were counted by the marker
  lane under its relatives, so it is listed beneath them (`counted_within`)
  with no share of its own and its own lane's reading kept as
  `secondary_percent`.
* **absorbed** — a rejected call's reads belong to the relative that took
  them in the competition, so the share the marker lane had read under the
  rejected name goes to that relative (`absorbed_percent`, `absorbed_from`
  on the relative; `reads_belong_to` on the rejected record). A relative
  with a marker reading of its own adds the share to it; one the marker
  catalogue had no entry for takes it as its share (basis `absorbed`), and
  a whole-genome estimate of the same population gives way to it. *Blautia
  massiliensis* at 5.5% in one sample was rejected because its reads mapped
  to *B. caecimuris* at higher identity; the 5.5% is *B. caecimuris*'s, and
  the composition adds up again.

`share_unification` in the inventory records the scale factor and the
counts. Nothing is summed across lanes: the primary lane's total is the
total, and the report's single list is sorted by this one share.

### Level against the reference

Two percentiles exist for an organism, and the report uses one of them.
`percentile` is prevalence-aware: a sample lacking the organism sits among
the reference people who also lack it, and a sample carrying it ranks
above every one of them, then by level among the carriers. That is the
right shape for the scoring machinery (an organism 9% of people carry is
unusual to have at all) and it is what GMWI and the FMT logic were
calibrated on, so it stays in `results.json`. It is the wrong number to
print beside a level: an organism carried by 9% of the reference at a
level a third of the typical carrier's read as "96th percentile, −61%".

`carrier_percentile` is the rank among reference people who carry the
organism, on the same readings whose median is `reference_percent`; the
report's tables, cards and groups print it and judge levels by it
(`Organism.level_percentile`). The deviation
(`reference_reading / reference_percent − 1`; +400% is five times the
typical carrier, −90% a tenth) is taken on that same reading and that same
median, so above the 50th percentile is always a positive deviation and
below always negative. The reading is the scoring lane's value rescaled to
its classified fraction - the number the rank was computed on - not the
composition share. Every rank names the lane whose reading it was taken on
(`reference_lane`, with that lane's raw reading as `reference_lane_percent`),
and where that reading and the share disagree fivefold the row says "not
comparable" instead of a rank. Rarity of carriage is stated on its own:
"carried by 9%", and an opportunist fewer than one in ten reference adults
carry is listed to watch whatever its level.

### What the report does and does not say

Detection statuses (`supported` / `provisional`), rejected calls with their
measurements, and the comparison with the previous reference set are in
`results.json` and in `results/EXPANSION_COMPARISON.md` (written beside the results by
`scripts/sample_comparison.py`). The report states
the organisms, their shares on one scale, their level against the
reference, how many methods saw each and the identifiers — not the
software's history.

Sibling GTDB species (`X` and `X_B`) placed by lanes that never both saw
the same organism — a marker bin labelled by the majority of its members,
a genome sketch placing the reads in the type cluster — fold into one
record shown as a complex (`Blautia_A wexlerae / wexlerae_B`), never as two
species. A lane that saw both siblings keeps them apart.

Status per organism:

* `supported` — two independent methods, or one method with a strong
  reading of its own, or competitive confirmation
* `provisional` — one method at a marginal reading; read classification
  alone is always provisional until confirmation
* `ambiguous` — a bin that maps to several R232 species (a complex)
* `rejected` — confirmation showed the reads belong to a relative; kept in
  `organism_inventory.rejected` with the reason

Counts are separate: named species, unnamed species-level clusters,
unresolved complexes, higher-rank findings; strains never enter species
richness. `organism_inventory.comparison` lists newly supported,
provisional, rejected and renamed organisms per sample.

### Percentiles per lane (`refs/expansion_cohort`)

A percentile is a position in a population measured the same way. The
scoring lane's population is the 3,027-adult curatedMetagenomicData cohort
(MetaPhlAn 3); an organism only the expanded lanes see cannot be placed in
it. `scripts/build_expansion_cohort.py` builds one population per lane by
running that lane's own engine over the same 100 adults of ENA study
PRJNA1271016 the pathway ranges rest on (3,000,000 read pairs per run,
prefix sampling; the manifest records depth, sampling, engine and database
release). `openbiota.expansion.cohorts` ranks a lane's readings against
its population with the scoring lane's arithmetic (prevalence-aware CLR
percentile). One population per organism: a lane rank is taken when the
scoring set never measured the organism, or when the scoring catalogue
reads it fivefold differently from the composition (its species boundary,
its markers) and the lane's own reading agrees with the share; every rank
carries `percentile_source` ("globdb cohort, n=100"). The catalogue marks
these ‡ and names the population. Populations are never merged or
averaged; exactly one is used per organism.

### Targeted assembly

The assembly stage runs when SingleM leaves lineages above species at
reconstructable coverage, and assembles only the reads Kraken2/UHGG left
unclassified or placed above species (`assembly.reads` records the
selection and its counts). Those reads are what an uncatalogued organism
is made of; assembling the 85% that already have a species would cost an
hour per sample and say nothing new. metaSPAdes runs `--only-assembler`
(BayesHammer asserted on this platform). Contigs are reported per lineage
with marker-family counts and are never counted as organisms.

### Strain analysis

StrainPhlAn 4 on the Jan26 consensus markers. Marker extraction is one
seqkit pass over the kept marker dump (`refs/metaphlan4_db/<index>.fna`)
using a marker→clade table derived once from the database pickle; the
output is byte-identical to `extract_markers.py`, which would re-dump the
21 GB index and parse it once per clade. Placed clades are written onto
their organism's `strain` record (nearest reference genome, distance,
markers); a placement never changes organism counts. Tool failures are
not cached, so the next run retries them.

## Operating points (initial; calibrated by `scripts/benchmark_lanes.py`)

* Kraken: `--confidence 0.15 --minimum-hit-groups 3`; inventory floor
  1,000 read pairs and 20,000 distinct minimizers; Bracken ≥ 10
* Confirmation: ≥ 100 uniquely mapping fragments, ≥ 3% breadth, ≥ 97%
  identity, evenness ≥ 0.25, with detected relatives competing. Evenness is
  observed breadth over the breadth the same mean depth would give if reads
  fell at random (1 − e^−depth): near 1 for reads spread over the genome,
  near 0 for reads piled into shared regions. (An earlier definition, 1 − CV
  of per-contig depth, scored every fragmented draft genome near 0 and kept
  organisms with tens of thousands of unique fragments provisional.)
* SingleM assembly trigger: unresolved lineage at ≥ 3× marker coverage
* StrainPhlAn: ≥ 20 reconstructed markers, every eligible organism the
  inventory lists (most abundant first), ≥ 3 reference genomes, a 4-hour
  budget per run for the trees (extraction and genome fetching excluded);
  each clade's result is cached, so a run that hits the budget records what
  it deferred and the next run resumes there
* Confirmation is relative as well as absolute: a same-genus relative that
  took more unique fragments at higher identity in the same competition
  rejects a candidate under 500 fragments ("the reads belong to the relative")
* Candidates: every expansion-only organism, every provisional or ambiguous
  call, and every single-method call with a same-genus relative detected -
  all of them, not a budgeted subset (the 400-candidate cap is a safety).
  Each lane supplies the genome its call rests on (the GlobDB genome behind
  a cluster, the panel genome a Kraken hit was classified to, a mOTUs
  cluster's representative, the UHGG species representative); a named GTDB
  species falls back to its R232 representative. Candidates with no genome
  anywhere (MetaPhlAn 3 names with no GTDB mapping, unresolved Jun23
  complexes) stay provisional and say so.

## Benchmarks (`scripts/benchmark_lanes.py`)

Truth communities are drawn from HRGM2 representatives frozen to their
R232 canonical species before simulation (InSilicoSeq, HiSeq model);
exact-reference and held-out-strain tiers are scored separately, a called
GlobDB/UHGG unit whose genome is the truth genome's species (skani ANI/AF)
is credited to that species, and parent or complex calls are never exact
true positives. Every lane is ablated against the installed baseline
(MetaPhlAn 3 + Jun23 + GTDB232 sylph), competitive confirmation runs inside
the loop, and precision, recall, false supported species and paired
recall gain carry Wilson or bootstrap 95% intervals.

`--negative N` adds organism-free backgrounds (human chr21/22, cattle and
chicken genome sequence) with the spec's gate of zero false supported
microbial species; provisional calls on those backgrounds are listed too.
Real blanks are not available and the summary says so.

## Completeness

The run's completeness gate has 24 required stages, including every lane,
confirmation, strain analysis, BacDive and assembly. A stage that did not
produce its result is printed as MISSING and the run exits 2.

## Licences

GlobDB CC BY-SA 4.0 · mOTUs data CC BY 4.0 · GTDB CC BY-SA 4.0 · BacDive
CC BY 4.0 (contact DSMZ for commercial use) · UHGG under EMBL-EBI terms ·
HRGM2 and ELGG CC BY 4.0 (Zenodo) · MetaPhlAn MIT, database per bioBakery.
