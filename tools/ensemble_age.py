"""Does ridge + TRPCA together beat ridge alone? Decided on held-out studies.

Two estimators with similar accuracy but different inductive biases often
blend well: ridge reads carriage of individual species linearly, TRPCA
attends over robust principal components. If their errors are less than
perfectly correlated, some weighted average beats either.

The weight is not allowed to see the test fold. For each held-out study the
blend weight is fitted on the out-of-fold predictions of the *other*
studies only, so the reported figure is what a new sample would get. A
weight chosen on the same predictions it is scored against would manufacture
an improvement out of nothing.

The ensemble ships only if it wins. That is the whole point of measuring.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from sweep_trpca import fit_once  # noqa: E402
from train_trpca import ADULT_MIN_AGE, SEED, metrics, ridge_fold  # noqa: E402

BEST = {"n_pcs": 48, "hidden": 128, "layers": 1, "lr": 1e-3}


def loso_predictions(X, y, g, fitter, label):
    pred = np.full(len(y), np.nan)
    t0 = time.monotonic()
    studies = sorted(set(g))
    for i, s in enumerate(studies, 1):
        held = np.where(g == s)[0]
        train = np.setdiff1d(np.arange(len(y)), held)
        if len(train) < 50:
            continue
        pred[held] = fitter(X[train], y[train], X[held])
        if i % 10 == 0 or i == len(studies):
            ok = np.isfinite(pred)
            print(f"  {label} [{i:2d}/{len(studies)}] {time.monotonic() - t0:5.0f}s "
                  f"MAE {np.mean(np.abs(pred[ok] - y[ok])):.2f}", flush=True)
    return pred


def best_weight(y, a, b, mask):
    """Weight w minimising MAE of w*a + (1-w)*b, on `mask` only."""
    grid = np.linspace(0.0, 1.0, 101)
    errs = [np.mean(np.abs((w * a[mask] + (1 - w) * b[mask]) - y[mask])) for w in grid]
    return float(grid[int(np.argmin(errs))])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", default="refs/cmd/age_training_matrix.npz")
    ap.add_argument("--all-ages", action="store_true")
    ap.add_argument("--out", default="refs/age/ENSEMBLE_V1")
    args = ap.parse_args()

    d = np.load(args.matrix, allow_pickle=True)
    X, age, study = d["X"], d["age"].astype(float), d["study"].astype(str)
    rows = np.ones(len(age), bool) if args.all_ages else age >= ADULT_MIN_AGE
    X, y, g = X[rows], age[rows], study[rows]
    print(f"{'all ages' if args.all_ages else 'adults'}: {len(y):,} across "
          f"{len(set(g))} studies\n", flush=True)

    ridge = loso_predictions(
        X, y, g, lambda A, ya, B: ridge_fold(A, ya, B)[0], "ridge"
    )
    trpca = loso_predictions(
        X, y, g,
        lambda A, ya, B: fit_once(
            A, ya, B, **BEST, wd=1e-4, epochs=200, batch=128, seed=SEED
        ),
        "trpca",
    )

    ok = np.isfinite(ridge) & np.isfinite(trpca)
    mr, mt = metrics(y[ok], ridge[ok]), metrics(y[ok], trpca[ok])
    corr = float(np.corrcoef(ridge[ok] - y[ok], trpca[ok] - y[ok])[0, 1])

    # Honest blend: the weight for each study comes from the other studies.
    blend = np.full(len(y), np.nan)
    weights = {}
    for s in sorted(set(g)):
        held = (g == s) & ok
        other = ok & ~(g == s)
        if not held.any() or not other.any():
            continue
        w = best_weight(y, ridge, trpca, other)
        weights[s] = w
        blend[held] = w * ridge[held] + (1 - w) * trpca[held]
    okb = np.isfinite(blend)
    mb = metrics(y[okb], blend[okb])

    print(f"\nridge    MAE {mr['mae_years']:.3f}  R2 {mr['r2']:+.3f}  slope {mr['calibration_slope']:.2f}")
    print(f"trpca    MAE {mt['mae_years']:.3f}  R2 {mt['r2']:+.3f}  slope {mt['calibration_slope']:.2f}")
    print(f"blend    MAE {mb['mae_years']:.3f}  R2 {mb['r2']:+.3f}  slope {mb['calibration_slope']:.2f}")
    print(f"\nresidual correlation between the two models: {corr:+.3f}")
    print(f"median blend weight on ridge: {np.median(list(weights.values())):.2f}")
    gain = mr["mae_years"] - mb["mae_years"]
    print(f"\nensemble improves MAE by {gain:+.3f} y "
          f"({'SHIP' if gain > 0.05 else 'DO NOT SHIP — no real gain'})")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out / "loso.npz", y=y, study=g, ridge=ridge, trpca=trpca, blend=blend)
    (out / "comparison.json").write_text(json.dumps({
        "scope": "all_ages" if args.all_ages else "adults",
        "folds": "leave-one-study-out; blend weight fitted on the other studies only",
        "ridge": mr, "trpca": mt, "blend": mb,
        "residual_correlation": corr,
        "median_ridge_weight": float(np.median(list(weights.values()))),
        "trpca_architecture": BEST,
        "verdict": "ship_ensemble" if gain > 0.05 else "keep_ridge_alone",
    }, indent=1))
    print(f"wrote {out}/comparison.json")


if __name__ == "__main__":
    main()
