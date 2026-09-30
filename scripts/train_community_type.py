#!/usr/bin/env python
"""Train and freeze the gut community-type model (BUILD_SPEC_v0.8.3 §6.2).

Run once. Writes `extension/community_type_model.json`, which the report
loads and never retrains. Nothing at report time clusters anything, and no
sample is ever ranked against the other samples in its own run.

    .venv/bin/python scripts/train_community_type.py

The procedure is the specification's, in order:

1. A frozen reference from one abundance method and a fixed genus universe,
   each sample normalised to one.
2. Jensen-Shannon divergence with natural logs; PAM on `sqrt(JSD)`.
3. `k = 2..6`, judged by study-stratified bootstrap stability and held-out
   silhouette. If no k is stable, the model is written with `stable: false`
   and assignment returns composition and distances rather than a category.
4. Medoids, features and the evaluation report are frozen together.
5. Groups are named for the genera that actually distinguish them.
"""

from __future__ import annotations

import argparse
import json
import random
import sys
import time
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from openbiota.extension import communitytype as CT  # noqa: E402

COHORT = REPO / "refs" / "taxonomic_cohort.json"

#: A k is only considered if at least this fraction of held-out samples keep
#: their assignment when the cohort is resampled by study.
STABILITY_THRESHOLD = 0.80

#: Bootstrap replicates for the stability estimate.
N_BOOTSTRAP = 20

#: Samples drawn for the distance matrix. PAM is O(n^2) in memory and the
#: cohort is 3,027 samples; a stratified subsample of this size keeps the
#: matrix tractable while preserving every study.
N_TRAIN = 1200

SEED = 20260923


def load_cohort() -> tuple[list[dict[str, float]], list[str]]:
    """Each reference sample as species -> percent, plus its study name."""
    raw = json.loads(COHORT.read_text(encoding="utf-8"))
    taxa: list[str] = raw["taxa"]
    matrix: list[list[float]] = raw["abundance"]  # taxa x samples
    sample_ids: list[str] = raw["sample_ids"]
    metadata: dict[str, dict[str, Any]] = raw["metadata"]
    samples: list[dict[str, float]] = []
    studies: list[str] = []
    for j, sample_id in enumerate(sample_ids):
        row = {taxa[i]: matrix[i][j] for i in range(len(taxa)) if matrix[i][j] > 0}
        if not row:
            continue
        samples.append(row)
        studies.append(str((metadata.get(sample_id) or {}).get("study_name") or "unknown"))
    return samples, studies


def stratified_subsample(studies: list[str], n: int, seed: int) -> list[int]:
    """Up to a proportional share from every study, so none is lost."""
    rng = random.Random(seed)
    by_study: dict[str, list[int]] = {}
    for i, study in enumerate(studies):
        by_study.setdefault(study, []).append(i)
    per_study = max(1, n // max(1, len(by_study)))
    chosen: list[int] = []
    for study in sorted(by_study):
        pool = by_study[study][:]
        rng.shuffle(pool)
        chosen.extend(pool[:per_study])
    rng.shuffle(chosen)
    return sorted(chosen[:n])


def distance_matrix(vectors: list[list[float]], *, chunk: int = 256) -> list[list[float]]:
    """Pairwise `sqrt(JSD)`, vectorised.

    The runtime module's scalar `root_jsd` is the definition; this is the
    same arithmetic in bulk, because training computes a hundred of these
    matrices and a pure-Python inner loop would take hours. The two are
    checked against each other in `verify_against_reference`.
    """
    import numpy as np  # noqa: PLC0415

    P = np.asarray(vectors, dtype=np.float64)
    n, _ = P.shape
    out = np.zeros((n, n), dtype=np.float64)
    # x log x, with 0 log 0 = 0, precomputed once per row.
    xlogx = np.zeros_like(P)
    positive = P > 0
    np.log(P, out=xlogx, where=positive)
    xlogx *= P
    row_entropy = xlogx.sum(axis=1)
    for start in range(0, n, chunk):
        stop = min(start + chunk, n)
        block = P[start:stop]                      # (b, d)
        m = 0.5 * (block[:, None, :] + P[None, :, :])   # (b, n, d)
        logm = np.zeros_like(m)
        np.log(m, out=logm, where=m > 0)
        cross = (m * logm).sum(axis=2)
        jsd = 0.5 * (row_entropy[start:stop, None] + row_entropy[None, :]) - cross
        out[start:stop] = np.sqrt(np.clip(jsd, 0.0, None))
    np.fill_diagonal(out, 0.0)
    out = 0.5 * (out + out.T)
    return out.tolist()


def verify_against_reference(vectors: list[list[float]], matrix: list[list[float]],
                             *, n_checks: int = 200, tolerance: float = 1e-9) -> None:
    """The bulk matrix must agree with the module's scalar definition."""
    rng = random.Random(SEED)
    n = len(vectors)
    worst = 0.0
    for _ in range(n_checks):
        i, j = rng.randrange(n), rng.randrange(n)
        worst = max(worst, abs(matrix[i][j] - CT.root_jsd(vectors[i], vectors[j])))
    if worst > tolerance:
        raise SystemExit(f"vectorised distance disagrees with the definition by {worst:.2e}")
    print(f"  distance matrix agrees with the scalar definition to {worst:.2e}", flush=True)


def bootstrap_stability(
    vectors: list[list[float]], studies: list[str], k: int, *, replicates: int, seed: int,
    base_distance: list[list[float]] | None = None,
) -> float:
    """How often a held-out sample keeps its assignment under study resampling.

    Studies are resampled with replacement, not samples: the question is
    whether the grouping survives a different mix of cohorts, which is the
    way this kind of model actually fails.
    """
    rng = random.Random(seed)
    unique_studies = sorted(set(studies))
    by_study: dict[str, list[int]] = {}
    for i, study in enumerate(studies):
        by_study.setdefault(study, []).append(i)

    base_d = base_distance if base_distance is not None else distance_matrix(vectors)
    base_medoids = CT.pam(base_d, k, seed=seed)
    base_labels = [
        min(range(k), key=lambda j: base_d[i][base_medoids[j]]) for i in range(len(vectors))
    ]
    agreements: list[float] = []
    for r in range(replicates):
        drawn = [rng.choice(unique_studies) for _ in unique_studies]
        idx = sorted({i for s in drawn for i in by_study[s]})
        if len(idx) < k * 3:
            continue
        sub_vectors = [vectors[i] for i in idx]
        sub_d = distance_matrix(sub_vectors)
        sub_medoids = CT.pam(sub_d, k, seed=seed + r + 1)
        # Held-out samples: those no replicate draw included.
        held_out = [i for i in range(len(vectors)) if i not in set(idx)]
        if not held_out:
            continue
        # Map each replicate medoid to the base medoid it is nearest to, so
        # the two labellings can be compared at all.
        mapping: dict[int, int] = {}
        for j, m in enumerate(sub_medoids):
            mv = sub_vectors[m]
            mapping[j] = min(range(k), key=lambda b: CT.root_jsd(mv, vectors[base_medoids[b]]))
        same = 0
        for i in held_out:
            sub_label = min(range(k), key=lambda j: CT.root_jsd(vectors[i], sub_vectors[sub_medoids[j]]))
            if mapping[sub_label] == base_labels[i]:
                same += 1
        agreements.append(same / len(held_out))
    return sum(agreements) / len(agreements) if agreements else 0.0


def held_out_silhouette(
    vectors: list[list[float]], k: int, *, seed: int, holdout: float = 0.30
) -> float:
    rng = random.Random(seed)
    idx = list(range(len(vectors)))
    rng.shuffle(idx)
    cut = int(len(idx) * (1 - holdout))
    train_idx, test_idx = sorted(idx[:cut]), sorted(idx[cut:])
    train_vectors = [vectors[i] for i in train_idx]
    medoids = CT.pam(distance_matrix(train_vectors), k, seed=seed)
    medoid_vectors = [train_vectors[m] for m in medoids]
    test_vectors = [vectors[i] for i in test_idx]
    labels = [
        min(range(k), key=lambda j: CT.root_jsd(v, medoid_vectors[j])) for v in test_vectors
    ]
    return CT.silhouette(distance_matrix(test_vectors), labels, k)


def name_groups(
    features: list[str], vectors: list[list[float]], labels: list[int], k: int
) -> tuple[list[str], list[tuple[str, ...]]]:
    """Name each group for the genera that actually distinguish it."""
    taxa = [CT.distinguishing(vectors, labels, features, c) for c in range(k)]
    names: list[str] = []
    used: set[str] = set()
    for c in range(k):
        top = taxa[c]
        base = f"{top[0]}-dominant" if top else f"group-{c + 1}"
        name = base
        n = 2
        while name in used:
            name = f"{base}-{n}"
            n += 1
        used.add(name)
        names.append(name)
    return names, taxa


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, default=CT.MODEL_FILE)
    parser.add_argument("--n-train", type=int, default=N_TRAIN)
    parser.add_argument("--bootstrap", type=int, default=N_BOOTSTRAP)
    args = parser.parse_args()

    started = time.monotonic()
    print(f"reading {COHORT}", flush=True)
    samples, studies = load_cohort()
    print(f"  {len(samples):,} reference samples, {len(set(studies))} studies", flush=True)

    keep = stratified_subsample(studies, args.n_train, SEED)
    train_samples = [samples[i] for i in keep]
    train_studies = [studies[i] for i in keep]
    print(f"  {len(train_samples):,} used for training, every study represented "
          f"({len(set(train_studies))})", flush=True)

    features, vectors = CT.genus_matrix(train_samples)
    print(f"  genus universe: {len(features)} genera", flush=True)

    full_distance = distance_matrix(vectors)
    verify_against_reference(vectors, full_distance)

    report = CT.TrainingReport()
    for k in range(2, 7):
        t = time.monotonic()
        sil = held_out_silhouette(vectors, k, seed=SEED)
        stab = bootstrap_stability(vectors, train_studies, k, replicates=args.bootstrap,
                                   seed=SEED, base_distance=full_distance)
        report.k_evaluated.append(k)
        report.silhouette[k] = round(sil, 6)
        report.stability[k] = round(stab, 6)
        print(f"  k={k}: held-out silhouette {sil:+.4f}, study-bootstrap stability {stab:.3f} "
              f"({time.monotonic() - t:.0f}s)", flush=True)

    eligible = [k for k in report.k_evaluated if report.stability[k] >= STABILITY_THRESHOLD]
    if eligible:
        report.selected_k = max(eligible, key=lambda k: report.silhouette[k])
        report.stable = True
        report.reason = (
            f"k={report.selected_k} had the highest held-out silhouette "
            f"({report.silhouette[report.selected_k]:+.4f}) among the k whose study-stratified "
            f"bootstrap stability reached {STABILITY_THRESHOLD:.0%}"
        )
    else:
        report.selected_k = max(report.k_evaluated, key=lambda k: report.stability[k])
        report.stable = False
        best = report.stability[report.selected_k]
        report.reason = (
            f"no k reached the {STABILITY_THRESHOLD:.0%} stability threshold (best was k="
            f"{report.selected_k} at {best:.3f}), so discrete types are not named; assignment "
            "returns the composition and its distances instead"
        )
    print(f"\n{report.reason}", flush=True)

    k = int(report.selected_k)
    full_d = full_distance
    medoid_idx = CT.pam(full_d, k, seed=SEED)
    labels = [min(range(k), key=lambda j: full_d[i][medoid_idx[j]]) for i in range(len(vectors))]
    names, taxa = name_groups(features, vectors, labels, k)

    # tau for soft similarities: the median nearest-medoid distance on
    # held-out data, so a typical sample's similarities are informative
    # rather than saturated at one.
    nearest = sorted(min(full_d[i][m] for m in medoid_idx) for i in range(len(vectors)))
    tau = max(1e-6, nearest[len(nearest) // 2])

    sizes = [sum(1 for lab in labels if lab == c) for c in range(k)]
    # How much of a reference sample the genus universe typically names. A
    # sample's own coverage means nothing without this: 40% is unremarkable
    # against a cohort median of 42% and alarming against one of 90%.
    coverages = sorted(
        sum(CT.composition(sample).shares.get(g, 0.0) for g in features)
        for sample in train_samples
    )
    mid = len(coverages) // 2
    reference_coverage = {
        "median": coverages[mid],
        "p10": coverages[int(len(coverages) * 0.10)],
        "p90": coverages[int(len(coverages) * 0.90)],
        "note": (
            "share of a reference sample's abundance that the genus universe names; a new "
            "sample's coverage is only interpretable beside this"
        ),
    }
    payload = {
        "schema_version": "openbiota.community-type-model/1.0",
        "spec_version": "0.8.3",
        "model_id": f"ext083.community-type/{k}-medoid/1.0",
        "built_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "method": {
            "divergence": "Jensen-Shannon, natural logs",
            "metric": "sqrt(JSD)",
            "clustering": "partitioning around medoids, deterministic seed",
            "feature_space": "genus shares over a fixed universe, renormalised",
            "seed": SEED,
        },
        "features": features,
        "k": k,
        "stable": report.stable,
        "tau": tau,
        "labels": names,
        "distinguishing_taxa": [list(t) for t in taxa],
        "medoids": [vectors[m] for m in medoid_idx],
        "group_sizes": sizes,
        "reference_feature_coverage": reference_coverage,
        "evaluation": report.to_json(),
        "cohort": {
            "source": "refs/taxonomic_cohort.json",
            "n_reference_samples": len(samples),
            "n_used_for_training": len(train_samples),
            "n_studies": len(set(studies)),
            "profiler": "MetaPhlAn 3.0, species level",
            "note": (
                "A frozen external cohort. A sample is never ranked against the other samples "
                "in its own run."
            ),
        },
        "limitations": [
            "Group membership is a resemblance to a reference composition, not a diagnosis.",
            "A Bacteroides-dominant result does not establish a high-protein or high-fat diet.",
            "Relative composition cannot establish absolute microbial load.",
        ],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=1), encoding="utf-8")
    print(f"\nwrote {args.out}", flush=True)
    print(f"  k={k} stable={report.stable} groups={names} sizes={sizes} tau={tau:.4f}", flush=True)
    print(f"  {time.monotonic() - started:.0f}s total", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
