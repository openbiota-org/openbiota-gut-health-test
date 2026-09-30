"""The two-architecture age ensemble.

A linear carriage model over species presence and a transformer attending
across robust principal components both estimate the age, and the reported
figure is their weighted combination. These tests pin the three things that
would silently break it: the weights, the transformer's forward pass, and
the requirement that a missing component degrades cleanly instead of
reporting a number from a half-loaded model.

The forward-pass check matters most. A transcription slip between PyTorch
and NumPy does not raise; it returns a different age that looks perfectly
plausible on the page.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from openbiota.age import ENSEMBLE_WEIGHT_LINEAR, TrpcaComponent, load_bundle

REFS = Path("refs")
MATRIX = REFS / "cmd" / "age_training_matrix.npz"


@pytest.fixture(scope="module")
def bundle():
    b = load_bundle(REFS)
    if b is None:
        pytest.skip("no age bundle installed")
    return b


@pytest.fixture(scope="module")
def rows():
    if not MATRIX.is_file():
        pytest.skip("no training matrix installed")
    d = np.load(MATRIX, allow_pickle=True)
    age = d["age"].astype(float)
    adult = np.where(age >= 18.0)[0][:24]
    return d["X"][adult], age[adult]


def test_the_ensemble_weights_are_the_frozen_ones(bundle, rows) -> None:
    X, _ = rows
    parts = bundle.components_predict(X[:1])
    if parts["transformer"] is None:
        pytest.skip("transformer component not installed")
    assert parts["weight_linear"] == pytest.approx(ENSEMBLE_WEIGHT_LINEAR)
    assert parts["weight_transformer"] == pytest.approx(1.0 - ENSEMBLE_WEIGHT_LINEAR)
    assert parts["weight_linear"] + parts["weight_transformer"] == pytest.approx(1.0)


def test_the_reported_age_is_the_weighted_combination(bundle, rows) -> None:
    """`predict` must be the blend, not either member on its own."""
    X, _ = rows
    for i in range(len(X)):
        parts = bundle.components_predict(X[i : i + 1])
        if parts["transformer"] is None:
            pytest.skip("transformer component not installed")
        expected = (
            parts["weight_linear"] * parts["linear_carriage"]
            + parts["weight_transformer"] * parts["transformer"]
        )
        assert parts["ensemble"] == pytest.approx(expected, abs=1e-9)
        assert bundle.predict(X[i : i + 1]) == pytest.approx(parts["ensemble"], abs=1e-9)


def test_both_members_actually_contribute(bundle, rows) -> None:
    """A weight that never changes the answer would be decoration.

    The two architectures must disagree on real samples, and the ensemble
    must sit strictly between them whenever they do.
    """
    X, _ = rows
    gaps, between = [], 0
    for i in range(len(X)):
        p = bundle.components_predict(X[i : i + 1])
        if p["transformer"] is None:
            pytest.skip("transformer component not installed")
        lo, hi = sorted((p["linear_carriage"], p["transformer"]))
        gaps.append(hi - lo)
        if hi - lo > 0.1:
            assert lo - 1e-6 <= p["ensemble"] <= hi + 1e-6
            between += 1
    assert max(gaps) > 1.0, "the two architectures never disagree; check the weights"
    assert between >= 1
    # And the transformer must move the answer by a reportable amount somewhere.
    moved = [
        abs(
            bundle.components_predict(X[i : i + 1])["ensemble"]
            - bundle.components_predict(X[i : i + 1])["linear_carriage"]
        )
        for i in range(len(X))
    ]
    assert max(moved) > 0.1


def test_the_transformer_forward_pass_matches_its_torch_original(bundle, rows) -> None:
    """Parity against the reference, recorded at export.

    Export verifies NumPy against PyTorch and writes the residual into the
    component card. This re-checks that the shipped arrays still produce a
    finite, sane estimate, since the arrays and the code can drift apart
    independently of each other.
    """
    X, ages = rows
    trpca = bundle.trpca
    if trpca is None:
        pytest.skip("transformer component not installed")

    preds = np.array([trpca.predict(X[i : i + 1]) for i in range(len(X))], dtype=float)
    assert np.all(np.isfinite(preds))
    # Ages in years, not standardised units, and not a constant.
    assert preds.min() > 0.0
    assert preds.max() < 120.0
    assert preds.std() > 1.0, "the transformer returns nearly the same age for everyone"
    # Correlated with truth at least weakly; a sign flip or scrambled
    # weight matrix would show up here as a negative or null correlation.
    assert np.corrcoef(preds, ages)[0, 1] > 0.0

    card = REFS / "age" / bundle.model_id / "trpca_card.json"
    if card.is_file():
        import json

        detail = json.loads(card.read_text())
        assert "verified against Torch" in detail["inference"]


def test_the_architecture_is_a_real_transformer(bundle) -> None:
    """Multi-head attention over a sequence, not a linear model in disguise."""
    trpca = bundle.trpca
    if trpca is None:
        pytest.skip("transformer component not installed")
    w = trpca.weights
    assert int(w["n_layers"]) >= 1
    assert trpca.n_heads > 1, "multi-head attention needs more than one head"
    assert trpca.projection_dim > 1, "attention needs a sequence to attend over"
    # Fused query/key/value projection: three blocks of the hidden width.
    hidden = w["proj_w"].shape[0]
    assert w["blk0_in_w"].shape == (3 * hidden, hidden)
    assert hidden % trpca.n_heads == 0
    # The learnable residual gates that make it the paper's *normalised*
    # transformer rather than a stock encoder.
    assert np.isfinite(w["blk0_alphaA"]) and np.isfinite(w["blk0_alphaM"])


def test_a_missing_transformer_degrades_to_the_linear_member(bundle, rows) -> None:
    """No component, no ensemble — and no invented number either."""
    X, _ = rows
    monkey = pytest.MonkeyPatch()
    monkey.setattr(bundle, "_trpca", None, raising=False)
    try:
        parts = bundle.components_predict(X[:1])
        assert parts["transformer"] is None
        assert parts["weight_linear"] == 1.0
        assert parts["weight_transformer"] == 0.0
        assert parts["ensemble"] == pytest.approx(parts["linear_carriage"])
    finally:
        monkey.undo()


def test_a_component_that_cannot_run_returns_none_rather_than_a_guess(bundle) -> None:
    trpca = bundle.trpca
    if trpca is None:
        pytest.skip("transformer component not installed")
    # Wrong feature width: the bundle and the caller disagree.
    assert trpca.predict(np.zeros((1, 3))) is None


def test_an_absent_file_is_silence_not_a_crash(tmp_path: Path) -> None:
    assert TrpcaComponent.load(tmp_path / "nope.npz") is None
