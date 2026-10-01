# Measured capability (generated)

Generated 2026-09-30T23:09-07:00 by `scripts/measured_capability_doc.py` from the newest benchmark summary and the per-sample comparison. Numbers here are measurements; the spec's targets are listed beside them and are not claimed where they were not met.

## Benchmark `overnight_v4` (2026-09-29T20:30)

2 simulated communities × 600 species × 8,000,000 read pairs (InSilicoSeq, HiSeq error model; truth = HRGM2 representatives frozen to GTDB R232).

| measure | installed baseline (MP3 + Jun23 + GTDB sylph) | all lanes | spec target |
|---|---|---|---|
| mean species recall | 58.5% | **94.2%** | ≥ 95% at 8M pairs |
| pooled precision (95% CI) | — | **92.6%** (91.0%–94.0%) | ≥ 99% |
| false supported species / community | — | 33.0 (1.5 sibling-split, 31.5 same-genus other) | ≤ 5 |
| precision at species-complex level (sibling splits credited) | — | 92.9% | — |
| paired recall gain (bootstrap 95% CI) | — | 35.7% (33.7% to 37.7%) | CI above zero |

Gates: precision_ge_0.99 = not met, recall_ge_0.95_at_8M = not met, false_supported_le_5 = not met, gain_ci_above_zero = met

Organism-free backgrounds (human chr21/22 + cattle + chicken DNA, 2 communities): **0 false supported microbial species**, 2 provisional (read classification alone, never supported by rule). Real blanks: none available; disclosed.

Communities run: 2 of the 20 the spec asks for; seeds and depth tiers not yet exhausted. The confidence intervals are correspondingly wide and the gates above are reported as measured, not as achieved at scale.

What the false supported species are: on every community measured so far, all of them sit in a genus that is present in the truth - a relative of a true species called beside it, most often by read classification and the marker lanes agreeing on a shared-sequence neighbour - and none is an unrelated organism. Competitive confirmation rejects the shadows it can test; widening it to every same-genus co-call, and calibrating the Kraken operating point per genus, is the next precision step and is not done.

## Per-sample change on the project's 6 specimens (unknown truth)

These are real stool samples with no known truth; the table shows what the expansion changed for each, not sensitivity or accuracy. Samples are identified by neutral codes.

| sample | organisms | supported | seen by installed baseline | newly supported | strains placed | rejected | provisional |
|---|---|---|---|---|---|---|---|
| SAMPLE1_A01 | 336 | 262 | 206 | 74 | 38 | 101 | 74 |
| SAMPLE2_A02 | 355 | 274 | 222 | 74 | 46 | 94 | 80 |
| SAMPLE3_A03 | 257 | 164 | 131 | 44 | 40 | 123 | 93 |
| SAMPLE4_A04 | 370 | 317 | 269 | 71 | 57 | 102 | 52 |
| SAMPLE5_A05 | 469 | 411 | 328 | 112 | 46 | 98 | 58 |
| SAMPLE6_A06 | 405 | 324 | 254 | 93 | 44 | 113 | 81 |

The per-organism lists behind this table are written to `results/EXPANSION_COMPARISON.md` by `scripts/sample_comparison.py` and stay with the results; they name organisms sample by sample and are not part of the documentation.
