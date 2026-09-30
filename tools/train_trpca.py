"""Train TRPCA (Myers et al. 2025) on the same folds as the shipped ridge model.

Transformer-based Robust Principal Component Analysis: RPCA vectors are
treated as a sequence, projected into several "views", and passed through a
normalised transformer encoder with multi-head attention. Architecture is
transcribed from the authors' reference implementation at
https://github.com/tydymy/TRPCA (`TRPCA/trpca.py`), which ships no trained
weights, so the model is trained here from our own cohort.

Two rules govern this script.

Identical folds. The comparison against ridge is worthless unless both see
exactly the same leave-one-study-out splits over exactly the same adult
subset, with every transformation — prevalence filter, RCLR, PCA, target
standardisation — fitted inside the fold. A single statistic computed on the
full matrix would leak the test study into the training set and flatter the
transformer.

Torch trains, NumPy serves. PyTorch has no wheel for this machine's Python,
so training runs in a separate interpreter and the fitted weights are
exported as plain arrays. The report's runtime keeps its single dependency,
and inference is a handful of matrix multiplications.

Run with the training interpreter, not the project one:

    .venv-trpca/bin/python tools/train_trpca.py
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import torch
from sklearn.decomposition import PCA
from torch import nn
from torch.nn import functional as F  # noqa: N812

# Matches openbiota.age
ADULT_MIN_AGE = 18.0
MIN_PREVALENCE = 0.02

# Reference defaults from the authors' example notebook.
NUM_PCS = 48
HIDDEN_DIM = 256
NUM_LAYERS = 1
PROJECTION_DIM = 4
NUM_HEADS = 4
BATCH_SIZE = 128
EPOCHS = 60
LR = 1e-3
WEIGHT_DECAY = 1e-4
SEED = 20260908


def rclr(X: np.ndarray) -> np.ndarray:
    """Robust centred log-ratio, identical to `openbiota.age.rclr`."""
    mask = X > 0
    with np.errstate(divide="ignore"):
        logs = np.where(mask, np.log(np.where(mask, X, 1.0)), 0.0)
    counts = mask.sum(axis=1, keepdims=True)
    counts = np.where(counts == 0, 1, counts)
    means = logs.sum(axis=1, keepdims=True) / counts
    return np.where(mask, logs - means, 0.0)


class NormalizedTransformerBlock(nn.Module):
    """Transcribed from the reference implementation.

    Every residual update is L2-normalised and scaled by a learnable gate,
    which is what the authors mean by "normalized transformer".
    """

    def __init__(self, dim: int, hidden: int) -> None:
        super().__init__()
        self.attention = nn.MultiheadAttention(dim, NUM_HEADS, dropout=0.0, batch_first=True)
        self.mlp = nn.Sequential(nn.Linear(dim, hidden), nn.ReLU(), nn.Linear(hidden, dim))
        self.alphaA = nn.Parameter(torch.tensor(1.0))
        self.alphaM = nn.Parameter(torch.tensor(1.0))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = F.normalize(x, p=2, dim=-1)
        hA, _ = self.attention(x, x, x)
        hA = F.normalize(hA, p=2, dim=-1)
        x = F.normalize(x + self.alphaA * (hA - x), p=2, dim=-1)
        hM = F.normalize(self.mlp(x), p=2, dim=-1)
        return F.normalize(x + self.alphaM * (hM - x), p=2, dim=-1)


def positional_encoding(length: int, dim: int) -> torch.Tensor:
    pe = torch.zeros(length, dim)
    pos = torch.arange(0, length, dtype=torch.float).unsqueeze(1)
    div = torch.exp(torch.arange(0, dim, 2).float() * (-np.log(10000.0) / dim))
    pe[:, 0::2] = torch.sin(pos * div)
    pe[:, 1::2] = torch.cos(pos * div)
    return pe.unsqueeze(0)


class TRPCA(nn.Module):
    def __init__(self, input_dim: int, hidden_dim: int, num_layers: int, projection_dim: int) -> None:
        super().__init__()
        self.projection_dim = projection_dim
        self.pca_projection = nn.Linear(input_dim, hidden_dim)
        self.view_generator = nn.Sequential(
            nn.Linear(hidden_dim, projection_dim * hidden_dim),
            nn.LayerNorm(projection_dim * hidden_dim),
        )
        self.register_buffer("pe", positional_encoding(projection_dim, hidden_dim))
        self.blocks = nn.ModuleList(
            NormalizedTransformerBlock(hidden_dim, hidden_dim * 2) for _ in range(num_layers)
        )
        self.regression_head = nn.Linear(hidden_dim, 1)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        b = x.shape[0]
        x = F.normalize(self.pca_projection(x), p=2, dim=-1)
        x = self.view_generator(x).view(b, self.projection_dim, -1)
        x = F.normalize(x, p=2, dim=-1) + self.pe
        for block in self.blocks:
            x = block(x)
        return self.regression_head(x.mean(dim=1)).squeeze(-1)


def fit_fold(
    Xtr: np.ndarray, ytr: np.ndarray, Xte: np.ndarray, *, seed: int, epochs: int = EPOCHS
) -> tuple[np.ndarray, dict]:
    """Fit every transformation and the network inside one fold."""
    torch.manual_seed(seed)
    np.random.seed(seed)

    keep = (Xtr > 0).mean(axis=0) >= MIN_PREVALENCE
    A, B = rclr(Xtr[:, keep]), rclr(Xte[:, keep])
    n_pcs = int(min(NUM_PCS, A.shape[0] - 1, A.shape[1]))
    pca = PCA(n_components=n_pcs, whiten=False, random_state=seed).fit(A)
    S, T = pca.transform(A), pca.transform(B)
    sd = S.std(axis=0) + 1e-9
    S, T = S / sd, T / sd

    # Standardise the target inside the fold too, so the loss scale does not
    # depend on the age distribution of whichever study is held out.
    ymu, ysd = float(ytr.mean()), float(ytr.std() + 1e-9)
    z = (ytr - ymu) / ysd

    model = TRPCA(n_pcs, HIDDEN_DIM, NUM_LAYERS, PROJECTION_DIM)
    opt = torch.optim.AdamW(model.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=epochs)
    Xt = torch.tensor(S, dtype=torch.float32)
    yt = torch.tensor(z, dtype=torch.float32)

    model.train()
    n = len(Xt)
    for _ in range(epochs):
        perm = torch.randperm(n)
        for i in range(0, n, BATCH_SIZE):
            idx = perm[i : i + BATCH_SIZE]
            opt.zero_grad()
            loss = F.smooth_l1_loss(model(Xt[idx]), yt[idx])
            loss.backward()
            opt.step()
        sched.step()

    model.eval()
    with torch.no_grad():
        pred = model(torch.tensor(T, dtype=torch.float32)).numpy() * ysd + ymu
    return pred, {"n_pcs": n_pcs, "n_features": int(keep.sum())}


def metrics(y: np.ndarray, pred: np.ndarray) -> dict:
    resid = pred - y
    ss_res = float(np.sum(resid**2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    slope = float(np.polyfit(pred, y, 1)[0]) if pred.std() > 1e-9 else 0.0
    return {
        "mae_years": float(np.mean(np.abs(resid))),
        "rmse_years": float(np.sqrt(np.mean(resid**2))),
        "r2": 1.0 - ss_res / ss_tot if ss_tot else float("nan"),
        "calibration_slope": slope,
        "n": int(len(y)),
    }


def ridge_fold(
    Xtr: np.ndarray, ytr: np.ndarray, Xte: np.ndarray, *, seed: int = 0, epochs: int = 0  # noqa: ARG001
) -> tuple[np.ndarray, dict]:
    """The shipped model: ridge on standardised species presence/absence.

    Reproduced here so both estimators are measured by one harness on one
    split, rather than compared across two scripts that differ in some
    detail nobody remembers.
    """
    keep = (Xtr > 0).mean(axis=0) >= MIN_PREVALENCE
    A = (Xtr[:, keep] > 0).astype(float)
    B = (Xte[:, keep] > 0).astype(float)
    mu, sd = A.mean(axis=0), A.std(axis=0) + 1e-9
    A, B = (A - mu) / sd, (B - mu) / sd
    n_pcs = int(min(200, A.shape[0] - 1, A.shape[1]))
    pca = PCA(n_components=n_pcs, random_state=0).fit(A)
    S, T = pca.transform(A), pca.transform(B)
    alpha = 1000.0
    G = S.T @ S + alpha * np.eye(S.shape[1])
    beta = np.linalg.solve(G, S.T @ (ytr - ytr.mean()))
    return T @ beta + ytr.mean(), {"n_pcs": n_pcs}


def evaluate(
    X: np.ndarray, y: np.ndarray, g: np.ndarray, *, fitter, split: str, epochs: int, seeds: int,
    label: str,
) -> tuple[np.ndarray, dict]:
    """Run one estimator under one split scheme."""
    pred = np.full(len(y), np.nan)
    t0 = time.monotonic()

    if split == "loso":
        # A new sample comes from a population the model has never seen.
        folds = [(s, np.where(g == s)[0]) for s in sorted(set(g))]
    else:
        # The paper's own scheme: a random 10% held out, stratified by study,
        # so every test sample's study is also in the training set.
        rng = np.random.default_rng(42)
        held: list[int] = []
        for s in sorted(set(g)):
            idx = np.where(g == s)[0]
            k = max(1, int(round(0.1 * len(idx))))
            held.extend(rng.choice(idx, size=k, replace=False).tolist())
        folds = [("stratified_holdout", np.array(sorted(held)))]

    for i, (name, hold) in enumerate(folds, 1):
        train = np.setdiff1d(np.arange(len(y)), hold)
        if len(train) < 50:
            continue
        runs = [
            fitter(X[train], y[train], X[hold], seed=SEED + k, epochs=epochs)[0]
            for k in range(seeds)
        ]
        pred[hold] = np.mean(runs, axis=0)
        if i % 5 == 0 or i == len(folds):
            ok = np.isfinite(pred)
            print(
                f"  {label}/{split} [{i:2d}/{len(folds)}] {time.monotonic() - t0:5.0f}s "
                f"running MAE {np.mean(np.abs(pred[ok] - y[ok])):.2f} y",
                flush=True,
            )
    ok = np.isfinite(pred)
    return pred, metrics(y[ok], pred[ok])


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", default="refs/cmd/age_training_matrix.npz")
    ap.add_argument("--out", default="refs/age/TRPCA_V1")
    ap.add_argument("--epochs", type=int, default=EPOCHS)
    ap.add_argument("--seeds", type=int, default=1, help="average this many seeds per fold")
    ap.add_argument("--split", default="loso", choices=["loso", "paper", "both"])
    ap.add_argument("--models", default="trpca", help="comma list: trpca,ridge")
    ap.add_argument("--all-ages", action="store_true",
                    help="include children; the paper's cohort spans 0-91 and the\n                          infant-versus-adult contrast is easy signal")
    args = ap.parse_args()

    d = np.load(args.matrix, allow_pickle=True)
    X, age, study = d["X"], d["age"].astype(float), d["study"].astype(str)
    keep_rows = np.ones(len(age), dtype=bool) if args.all_ages else (age >= ADULT_MIN_AGE)
    X, y, g = X[keep_rows], age[keep_rows], study[keep_rows]
    who = "all ages" if args.all_ages else "adults"
    print(f"{who} {len(y):,} across {len(set(g))} studies, {X.shape[1]} taxa "
          f"(age {y.min():.1f}-{y.max():.1f})", flush=True)

    fitters = {"trpca": fit_fold, "ridge": ridge_fold}
    splits = ["loso", "paper"] if args.split == "both" else [args.split]
    results: dict[str, dict] = {}
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    for model in [m.strip() for m in args.models.split(",") if m.strip()]:
        for split in splits:
            pred, m = evaluate(
                X, y, g, fitter=fitters[model], split=split,
                epochs=args.epochs, seeds=args.seeds, label=model,
            )
            key = f"{model}_{split}"
            results[key] = m
            np.savez_compressed(out / f"pred_{key}.npz", pred=pred, y=y, study=g)
            print(f"\n{key}: MAE {m['mae_years']:.2f} y | R2 {m['r2']:+.3f} | "
                  f"slope {m['calibration_slope']:.2f} | n {m['n']}\n", flush=True)

    (out / "loso_metrics.json").write_text(
        json.dumps(
            {
                "model": "TRPCA (Myers et al. 2025), transcribed from tydymy/TRPCA",
                "architecture": {
                    "num_pcs": NUM_PCS, "hidden_dim": HIDDEN_DIM,
                    "num_layers": NUM_LAYERS, "projection_dim": PROJECTION_DIM,
                    "num_heads": NUM_HEADS, "epochs": args.epochs, "seeds": args.seeds,
                },
                "splits": {
                    "loso": "leave-one-study-out; the deployment question",
                    "paper": "random 10% stratified by study, the authors' own scheme; "
                             "train and test share studies",
                },
                "results": results,
            },
            indent=1,
        )
    )
    print(f"wrote {out}/loso_metrics.json")


if __name__ == "__main__":
    main()
