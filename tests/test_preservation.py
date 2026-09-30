"""BUILD_SPEC_v0.8.3 section 1: the extension may not move an existing number.

AT001-AT006. The promise is that every existing score, label, component and
calculation fingerprint is identical for identical resolved inputs, with the
extension present or absent. These tests hold the current application against
frozen baselines captured before any v0.8.3 work, so a regression in a
protected module fails here rather than in a rendered report.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from pathlib import Path
from typing import Any

import pytest

from openbiota.extension import preservation as P

REPO = Path(__file__).resolve().parents[1]
FIXTURES = REPO / "tests" / "fixtures" / "preservation"
RESULTS = REPO / "results"


def _manifest() -> dict[str, Any]:
    path = FIXTURES / "protected_manifest.json"
    if not path.is_file():
        pytest.skip("preservation fixtures are not present in this checkout")
    return json.loads(path.read_text(encoding="utf-8"))


def _baseline(sample: str) -> dict[str, Any]:
    path = FIXTURES / f"{sample}.protected.json.gz"
    if not path.is_file():
        pytest.skip(f"no frozen baseline for {sample}")
    return json.loads(gzip.decompress(path.read_bytes()).decode("utf-8"))


def _current(sample: str) -> dict[str, Any]:
    from openbiota.samples import local_name, results_file

    path = results_file(sample)
    if not path.is_file():
        pytest.skip(f"{sample} has not been run in this checkout")
    return P.anonymise(json.loads(path.read_text(encoding="utf-8")), local=local_name(sample), published=sample)


SAMPLES = sorted(s for s in (_manifest()["samples"] if (FIXTURES / "protected_manifest.json").is_file() else {}))


@pytest.mark.parametrize("sample", SAMPLES)
def test_at001_at003_protected_results_equal_the_frozen_baseline(sample: str) -> None:
    """AT001/AT002/AT003: every protected value is unchanged.

    Disease scores, GMWI2, donor matching, strain analysis, pathogen
    screening, biofilm, mycobiome and the functional panels all live under the
    protected keys, so one comparison covers them.
    """
    report = P.compare(_baseline(sample), _current(sample))
    assert report.compared_leaves > 10_000, report.summary()
    assert report.preserved, (
        report.summary()
        + "\n"
        + "\n".join(
            f"  {d.path}: {d.baseline!r} -> {d.candidate!r} ({d.reason})"
            for d in report.differences[:25]
        )
    )


@pytest.mark.parametrize("sample", SAMPLES)
def test_at003_every_captured_metric_id_still_exists(sample: str) -> None:
    """A protected metric may not quietly disappear.

    Comparing values catches a changed number. It does not catch a whole panel
    being dropped, because a missing key has no value to compare - so the
    manifest's ID list is checked separately.
    """
    expected = _manifest()["samples"][sample]["metric_ids"]
    actual = P.metric_ids(_current(sample))
    for group, ids in expected.items():
        assert group in actual, f"{sample}: the whole {group} group is gone"
        missing = sorted(set(ids) - set(actual[group]))
        assert not missing, f"{sample}: {group} lost {len(missing)} ID(s): {missing[:8]}"


def test_at004_a_new_measurement_cannot_enter_a_protected_object() -> None:
    """AT004: additions live outside the protected keys.

    The check that makes this enforceable rather than aspirational: the
    extension writes one new top-level object, and `PROTECTED_KEYS` does not
    contain it. A future module that wrote into `panels` would fail the
    comparison above; this test states the intent for the reader.
    """
    assert "extension" not in P.PROTECTED_KEYS
    assert "ext083" not in P.PROTECTED_KEYS
    for sample in SAMPLES[:1]:
        current = _current(sample)
        for key in P.PROTECTED_KEYS:
            if key in current:
                blob = json.dumps(current[key])
                assert "openbiota.extension-metric/" not in blob, (
                    f"an extension metric was written inside the protected object {key!r}"
                )


@pytest.mark.parametrize("sample", SAMPLES)
def test_frozen_fixture_matches_its_recorded_hash(sample: str) -> None:
    """A fixture that silently changed would make the whole suite meaningless."""
    recorded = _manifest()["samples"][sample].get("frozen_protected_sha256")
    if not recorded:
        pytest.skip("no hash recorded for this fixture")
    raw = gzip.decompress((FIXTURES / f"{sample}.protected.json.gz").read_bytes())
    assert hashlib.sha256(raw).hexdigest() == recorded


def test_the_comparison_actually_detects_a_changed_value() -> None:
    """A preservation check that cannot fail proves nothing.

    Every protected type is tampered with in turn - a number, a label, a list
    length and a removed key - and each must be reported.
    """
    baseline = {
        "panels": [{"name": "but", "copies_per_100_genomes": 12.5, "confidence": "high"}],
        "biofilm": {"cards": [{"card": "H-C", "raw_index": 0.42}]},
        "standing_caveats": ["one", "two"],
    }
    assert P.compare(baseline, baseline).preserved

    changed_number = json.loads(json.dumps(baseline))
    changed_number["panels"][0]["copies_per_100_genomes"] = 12.6
    assert not P.compare(baseline, changed_number).preserved

    changed_label = json.loads(json.dumps(baseline))
    changed_label["panels"][0]["confidence"] = "low"
    assert not P.compare(baseline, changed_label).preserved

    shorter_list = json.loads(json.dumps(baseline))
    shorter_list["standing_caveats"] = ["one"]
    assert not P.compare(baseline, shorter_list).preserved

    removed_key = json.loads(json.dumps(baseline))
    del removed_key["biofilm"]["cards"][0]["raw_index"]
    report = P.compare(baseline, removed_key)
    assert not report.preserved
    assert report.differences[0].reason == "key removed"

    dropped_object = json.loads(json.dumps(baseline))
    del dropped_object["biofilm"]
    report = P.compare(baseline, dropped_object)
    assert not report.preserved
    assert report.protected_keys_missing == ("biofilm",)


def test_float_noise_is_tolerated_but_a_real_change_is_not() -> None:
    """Last-bit float drift is not a changed calculation; 0.1% is."""
    baseline = {"panels": [{"name": "but", "copies_per_100_genomes": 1.0}]}
    noisy = {"panels": [{"name": "but", "copies_per_100_genomes": 1.0 + 1e-13}]}
    real = {"panels": [{"name": "but", "copies_per_100_genomes": 1.001}]}
    assert P.compare(baseline, noisy).preserved
    assert not P.compare(baseline, real).preserved


def test_volatile_run_metadata_is_enumerated_not_pattern_matched() -> None:
    """Timings differ between runs; a measurement that mentions time does not.

    The exclusion list is exact paths. A rule like "skip anything containing
    'time'" would also excuse a changed value whose label mentions timing,
    which is exactly the kind of hole this contract cannot have.
    """
    assert P._is_volatile("run.total wall clock")
    assert P._is_volatile("run.stage timings.pathogen screen")
    assert not P._is_volatile("panels[0].elapsed_s")
    assert not P._is_volatile("biofilm.cards[0].raw_index")
    # Nothing under a protected scientific object may be excused.
    allowed_roots = {
        "run", "input", "diagnostics", "search", "pathogens", "mycobiome", "biofilm",
        "extended_catalogue", "genome_profile", "profile_similarity", "sequencing_quality",
    }
    for path in P.VOLATILE_PATHS:
        assert path.split(".")[0].split("[")[0] in allowed_roots, path
    # Inside a protected object, only a timing may be excused. Every such
    # entry ends in a timing field, so no measurement can hide behind one.
    for path in P.VOLATILE_PATHS:
        if path.split(".")[0] not in {"run", "input", "diagnostics", "search"}:
            assert path.endswith(("timings_s", "elapsed_s")), (
                f"{path} excuses more than a stage timing"
            )
