"""Train the production TRPCA component and export it as plain arrays.

The transformer trains here, in an interpreter that has PyTorch, and is then
frozen to a `.npz` of weights plus the preprocessing it was fitted with:
prevalence mask, PCA basis, component scaling, and the target
standardisation. `openbiota.age` reloads those arrays and runs the forward pass
in NumPy, so the report's runtime keeps its single dependency and inference
is a handful of matrix multiplications.

The exported arrays are checked against the Torch model on real training
rows before anything is written. A silent transcription error in the
forward pass would not crash; it would quietly return a different age.

    .venv-trpca/bin/python tools/export_trpca.py
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import torch
from sklearn.decomposition import PCA
from torch import nn
from torch.nn import functional as F  # noqa: N812

sys.path.insert(0, str(Path(__file__).parent))
from train_trpca import (  # noqa: E402
    ADULT_MIN_AGE,
    MIN_PREVALENCE,
    PROJECTION_DIM,
    SEED,
    TRPCA,
    rclr,
)

BEST = {"n_pcs": 48, "hidden": 128, "layers": 1, "lr": 1e-3}
EPOCHS = 200
BATCH = 128
WD = 1e-4
VAL_FRAC = 0.2
PATIENCE = 10


def train_full(X: np.ndarray, y: np.ndarray, seed: int) -> tuple[nn.Module, dict]:
    """Fit preprocessing and network on the whole adult training set."""
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)

    keep = (X > 0).mean(axis=0) >= MIN_PREVALENCE
    A = rclr(X[:, keep])
    k = int(min(BEST["n_pcs"], A.shape[0] - 1, A.shape[1]))
    pca = PCA(n_components=k, random_state=seed).fit(A)
    S = pca.transform(A)
    sd = S.std(axis=0) + 1e-9
    S = S / sd
    ymu, ysd = float(y.mean()), float(y.std() + 1e-9)
    z = (y - ymu) / ysd

    n = len(S)
    perm = rng.permutation(n)
    n_val = max(1, int(VAL_FRAC * n))
    va, tr = perm[:n_val], perm[n_val:]
    Xt = torch.tensor(S[tr], dtype=torch.float32)
    yt = torch.tensor(z[tr], dtype=torch.float32)
    Xv = torch.tensor(S[va], dtype=torch.float32)
    yv = torch.tensor(z[va], dtype=torch.float32)

    model = TRPCA(k, BEST["hidden"], BEST["layers"], PROJECTION_DIM)
    opt = torch.optim.AdamW(model.parameters(), lr=BEST["lr"], weight_decay=WD)
    best, state, bad, best_epoch = float("inf"), None, 0, 0
    for epoch in range(EPOCHS):
        model.train()
        order = torch.randperm(len(Xt))
        for i in range(0, len(Xt), BATCH):
            idx = order[i : i + BATCH]
            opt.zero_grad()
            F.smooth_l1_loss(model(Xt[idx]), yt[idx]).backward()
            opt.step()
        model.eval()
        with torch.no_grad():
            v = float(F.l1_loss(model(Xv), yv))
        if v < best - 1e-4:
            best, bad, best_epoch = v, 0, epoch
            state = {k2: t.detach().clone() for k2, t in model.state_dict().items()}
        else:
            bad += 1
            if bad >= PATIENCE:
                break
    if state is not None:
        model.load_state_dict(state)
    model.eval()
    return model, {
        "keep": keep, "pca_mean": pca.mean_, "pca_components": pca.components_,
        "score_sd": sd, "y_mean": ymu, "y_sd": ysd, "n_pcs": k,
        "val_l1_standardised": best, "best_epoch": best_epoch,
    }


def numpy_forward(w: dict, S: np.ndarray) -> np.ndarray:
    """The forward pass `openbiota.age` will run. Kept here to verify parity."""
    def l2(x, axis=-1):
        return x / (np.linalg.norm(x, axis=axis, keepdims=True) + 1e-12)

    x = l2(S @ w["proj_w"].T + w["proj_b"])
    x = x @ w["view_w"].T + w["view_b"]
    # LayerNorm over the flattened view vector
    mu = x.mean(axis=-1, keepdims=True)
    sd = x.std(axis=-1, keepdims=True)
    x = (x - mu) / (sd + 1e-5) * w["view_ln_w"] + w["view_ln_b"]
    b = x.shape[0]
    x = l2(x.reshape(b, PROJECTION_DIM, -1)) + w["pe"]

    n_layers = int(w["n_layers"])
    heads = int(w["n_heads"])
    for i in range(n_layers):
        p = f"blk{i}_"
        x = l2(x)
        d = x.shape[-1]
        qkv = x @ w[p + "in_w"].T + w[p + "in_b"]          # [b, L, 3d]
        q, k, v = np.split(qkv, 3, axis=-1)
        hd = d // heads
        def split(t):
            return t.reshape(b, -1, heads, hd).transpose(0, 2, 1, 3)
        q, k, v = split(q), split(k), split(v)
        att = q @ k.transpose(0, 1, 3, 2) / np.sqrt(hd)
        att = att - att.max(axis=-1, keepdims=True)
        att = np.exp(att)
        att = att / att.sum(axis=-1, keepdims=True)
        h = (att @ v).transpose(0, 2, 1, 3).reshape(b, -1, d)
        h = h @ w[p + "out_w"].T + w[p + "out_b"]
        h = l2(h)
        x = l2(x + w[p + "alphaA"] * (h - x))
        m = np.maximum(x @ w[p + "fc1_w"].T + w[p + "fc1_b"], 0.0)
        m = m @ w[p + "fc2_w"].T + w[p + "fc2_b"]
        m = l2(m)
        x = l2(x + w[p + "alphaM"] * (m - x))

    pooled = x.mean(axis=1)
    return (pooled @ w["head_w"].T + w["head_b"]).squeeze(-1)


def extract(model: nn.Module) -> dict:
    sd = model.state_dict()
    w: dict = {
        "proj_w": sd["pca_projection.weight"].numpy(),
        "proj_b": sd["pca_projection.bias"].numpy(),
        "view_w": sd["view_generator.0.weight"].numpy(),
        "view_b": sd["view_generator.0.bias"].numpy(),
        "view_ln_w": sd["view_generator.1.weight"].numpy(),
        "view_ln_b": sd["view_generator.1.bias"].numpy(),
        "pe": sd["pe"].numpy(),
        "head_w": sd["regression_head.weight"].numpy(),
        "head_b": sd["regression_head.bias"].numpy(),
        "n_layers": np.array(len(model.blocks)),
        "n_heads": np.array(4),
    }
    for i, _ in enumerate(model.blocks):
        p = f"blocks.{i}."
        q = f"blk{i}_"
        w[q + "in_w"] = sd[p + "attention.in_proj_weight"].numpy()
        w[q + "in_b"] = sd[p + "attention.in_proj_bias"].numpy()
        w[q + "out_w"] = sd[p + "attention.out_proj.weight"].numpy()
        w[q + "out_b"] = sd[p + "attention.out_proj.bias"].numpy()
        w[q + "fc1_w"] = sd[p + "mlp.0.weight"].numpy()
        w[q + "fc1_b"] = sd[p + "mlp.0.bias"].numpy()
        w[q + "fc2_w"] = sd[p + "mlp.2.weight"].numpy()
        w[q + "fc2_b"] = sd[p + "mlp.2.bias"].numpy()
        w[q + "alphaA"] = sd[p + "alphaA"].numpy()
        w[q + "alphaM"] = sd[p + "alphaM"].numpy()
    return w


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--matrix", default="refs/cmd/age_training_matrix.npz")
    ap.add_argument("--out", default="refs/age/MPA_WGS_AGE_V1/trpca.npz")
    args = ap.parse_args()

    d = np.load(args.matrix, allow_pickle=True)
    X, age, taxa = d["X"], d["age"].astype(float), d["taxa"].astype(str)
    adult = age >= ADULT_MIN_AGE
    X, y = X[adult], age[adult]
    print(f"training on {len(y):,} adults, {X.shape[1]} taxa", flush=True)

    model, pre = train_full(X, y, SEED)
    w = extract(model)

    # ---- parity: NumPy must match Torch on real rows --------------------- #
    A = rclr(X[:, pre["keep"]])
    S = ((A - pre["pca_mean"]) @ pre["pca_components"].T) / pre["score_sd"]
    with torch.no_grad():
        ref = model(torch.tensor(S[:256], dtype=torch.float32)).numpy()
    got = numpy_forward(w, S[:256].astype(np.float64))
    err = float(np.max(np.abs(ref - got)))
    print(f"numpy/torch parity: max |diff| = {err:.2e} (standardised units)")
    if err > 1e-3:
        raise SystemExit(f"forward-pass transcription mismatch: {err:.3e}")

    pred = numpy_forward(w, S.astype(np.float64)) * pre["y_sd"] + pre["y_mean"]
    print(f"in-sample MAE {np.mean(np.abs(pred - y)):.2f} y "
          f"(in-sample, not a generalisation estimate)")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        **w,
        "keep": pre["keep"],
        "taxa": taxa,
        "pca_mean": pre["pca_mean"],
        "pca_components": pre["pca_components"],
        "score_sd": pre["score_sd"],
        "y_mean": np.array(pre["y_mean"]),
        "y_sd": np.array(pre["y_sd"]),
        "projection_dim": np.array(PROJECTION_DIM),
    }
    np.savez_compressed(out, **payload)
    (out.parent / "trpca_card.json").write_text(json.dumps({
        "component": "TRPCA — transformer-based robust PCA age estimator",
        "source": "Myers et al., Communications Biology 8:1159 (2025); "
                  "architecture from github.com/tydymy/TRPCA",
        "architecture": {**BEST, "projection_dim": PROJECTION_DIM, "heads": 4},
        "training_n": int(len(y)),
        "n_pcs": int(pre["n_pcs"]),
        "best_epoch": int(pre["best_epoch"]),
        "inference": "NumPy forward pass; verified against Torch to "
                     f"{err:.1e} in standardised units",
    }, indent=1))
    print(f"wrote {out} ({out.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
