"""Give TRPCA its best shot before concluding anything about it.

The reference implementation tunes learning rate, weight decay and optimiser
with Optuna against an inner validation split, and defaults to 20 epochs
rather than the 60 a first guess suggests. Comparing a tuned ridge against an
untuned transformer would be a rigged fight, so this sweeps the transformer's
hyperparameters on the authors' own split scheme and reports the best.

Selection uses an inner validation split carved out of the training data
only. The held-out test fold is never consulted while choosing.
"""

from __future__ import annotations

import argparse
import itertools
import json
import sys
import time
from pathlib import Path

import numpy as np
import torch
from sklearn.decomposition import PCA
from torch.nn import functional as F  # noqa: N812

sys.path.insert(0, str(Path(__file__).parent))
from train_trpca import (  # noqa: E402
    ADULT_MIN_AGE,
    MIN_PREVALENCE,
    PROJECTION_DIM,
    SEED,
    TRPCA,
    metrics,
    rclr,
    ridge_fold,
)


def fit_once(
    Xtr, ytr, Xte, *, n_pcs, hidden, layers, lr, wd, epochs, batch, seed, val_frac=0.2,
):
    """One fit with early stopping on an inner validation split."""
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)

    keep = (Xtr > 0).mean(axis=0) >= MIN_PREVALENCE
    A, B = rclr(Xtr[:, keep]), rclr(Xte[:, keep])
    k = int(min(n_pcs, A.shape[0] - 1, A.shape[1]))
    pca = PCA(n_components=k, random_state=seed).fit(A)
    S, T = pca.transform(A), pca.transform(B)
    sd = S.std(axis=0) + 1e-9
    S, T = S / sd, T / sd
    ymu, ysd = float(ytr.mean()), float(ytr.std() + 1e-9)
    z = (ytr - ymu) / ysd

    n = len(S)
    perm = rng.permutation(n)
    n_val = max(1, int(val_frac * n))
    va, tr = perm[:n_val], perm[n_val:]
    Xt = torch.tensor(S[tr], dtype=torch.float32)
    yt = torch.tensor(z[tr], dtype=torch.float32)
    Xv = torch.tensor(S[va], dtype=torch.float32)
    yv = torch.tensor(z[va], dtype=torch.float32)

    model = TRPCA(k, hidden, layers, PROJECTION_DIM)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=wd)
    best, best_state, patience = float("inf"), None, 0
    for _ in range(epochs):
        model.train()
        order = torch.randperm(len(Xt))
        for i in range(0, len(Xt), batch):
            idx = order[i : i + batch]
            opt.zero_grad()
            F.smooth_l1_loss(model(Xt[idx]), yt[idx]).backward()
            opt.step()
        model.eval()
        with torch.no_grad():
            v = float(F.l1_loss(model(Xv), yv))
        if v < best - 1e-4:
            best, patience = v, 0
            best_state = {k2: t.detach().clone() for k2, t in model.state_dict().items()}
        else:
            patience += 1
            if patience >= 10:
                break
    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        return model(torch.tensor(T, dtype=torch.float32)).numpy() * ysd + ymu


def paper_split(g, seed=42):
    rng = np.random.default_rng(seed)
    held = []
    for s in sorted(set(g)):
        idx = np.where(g == s)[0]
        n = max(1, int(round(0.1 * len(idx))))
        held.extend(rng.choice(idx, size=n, replace=False).tolist())
    return np.array(sorted(held))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", default="refs/cmd/age_training_matrix.npz")
    ap.add_argument("--all-ages", action="store_true")
    args = ap.parse_args()

    d = np.load(args.matrix, allow_pickle=True)
    X, age, study = d["X"], d["age"].astype(float), d["study"].astype(str)
    rows = np.ones(len(age), bool) if args.all_ages else age >= ADULT_MIN_AGE
    X, y, g = X[rows], age[rows], study[rows]
    hold = paper_split(g)
    train = np.setdiff1d(np.arange(len(y)), hold)
    print(f"{'all ages' if args.all_ages else 'adults'}: {len(y):,} samples, "
          f"{len(hold)} held out\n", flush=True)

    base = ridge_fold(X[train], y[train], X[hold])[0]
    bm = metrics(y[hold], base)
    print(f"ridge baseline: MAE {bm['mae_years']:.2f} | R2 {bm['r2']:+.3f} | "
          f"slope {bm['calibration_slope']:.2f}\n", flush=True)

    grid = list(itertools.product(
        [48, 128, 256],      # n_pcs
        [128, 256],          # hidden
        [1, 2],              # layers
        [1e-3, 3e-4],        # lr
    ))
    results = []
    t0 = time.monotonic()
    for i, (n_pcs, hidden, layers, lr) in enumerate(grid, 1):
        pred = fit_once(
            X[train], y[train], X[hold], n_pcs=n_pcs, hidden=hidden, layers=layers,
            lr=lr, wd=1e-4, epochs=200, batch=128, seed=SEED,
        )
        m = metrics(y[hold], pred)
        results.append({"n_pcs": n_pcs, "hidden": hidden, "layers": layers, "lr": lr, **m})
        print(f"  [{i:2d}/{len(grid)}] pcs={n_pcs:3d} h={hidden:3d} L={layers} "
              f"lr={lr:.0e} -> MAE {m['mae_years']:5.2f} R2 {m['r2']:+.3f} "
              f"slope {m['calibration_slope']:.2f}  ({time.monotonic() - t0:.0f}s)", flush=True)

    results.sort(key=lambda r: r["mae_years"])
    best = results[0]
    print(f"\nbest TRPCA: {json.dumps({k: best[k] for k in ('n_pcs','hidden','layers','lr','mae_years','r2','calibration_slope')}, indent=1)}")
    print(f"ridge:      MAE {bm['mae_years']:.2f} | R2 {bm['r2']:+.3f}")
    Path("/tmp/trpca_sweep.json").write_text(
        json.dumps({"ridge": bm, "trpca": results}, indent=1)
    )


if __name__ == "__main__":
    main()
