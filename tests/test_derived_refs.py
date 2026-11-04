"""The bundles of reference data that travel with the repository.

Two kinds of test. The first kind looks at the committed bundles themselves:
each matches its manifest, holds what the module says it holds, and every
part is small enough for GitHub to take as an ordinary file. The second
kind exercises pack / restore / verify / check on a throwaway repository
under tmp_path, including the split-and-join path.
"""

from __future__ import annotations

import importlib.util
import io
import json
import os
import shutil
import sys
import tarfile
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts" / "derived_refs.py"
GITHUB_WARNING = 50 * 1024 * 1024
GITHUB_LIMIT = 100 * 1024 * 1024


def _load(module_name: str = "derived_refs"):
    spec = importlib.util.spec_from_file_location(module_name, SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    # Dataclasses resolve string annotations through sys.modules, so a
    # module loaded by path has to be registered before it executes.
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def fake_repo(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    """The module pointed at a throwaway repository with small fake refs."""
    module = _load("derived_refs_under_test")
    root = tmp_path / "repo"
    (root / "refs" / "age").mkdir(parents=True)
    (root / "refs" / "crosswalk").mkdir()
    (root / "refs" / "panels" / "mechanisms").mkdir(parents=True)
    (root / "refs" / "age" / "model.json").write_text('{"model": 1}')
    (root / "refs" / "crosswalk" / "gtdb.summary.json").write_text("{}")
    (root / "refs" / "panels" / "mechanisms" / "mech.bft1.fasta").write_text(">x\nACGT\n")
    (root / "refs" / "reference_ranges.json").write_text('{"ranges": []}')
    monkeypatch.setattr(module, "REPO", root)
    monkeypatch.setattr(module, "ARCHIVE_DIR", root / "data")
    items = (
        module.Item("refs/reference_ranges.json", "ranges"),
        module.Item("refs/age", "age model"),
        module.Item("refs/crosswalk", "crosswalks"),
        module.Item("refs/panels", "panels"),
        module.Item("refs/bacdive", "cache", required=False),
    )
    bundle = module.Bundle("derived-refs", lambda: items, "test bundle")
    monkeypatch.setattr(module, "BUNDLES", {"derived-refs": bundle})
    monkeypatch.setattr(module, "DEFAULT_ORDER", ("derived-refs",))
    return module, root, bundle


# --------------------------------------------------------------------------- #
# the committed bundles
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("name", ["derived-refs", "uniprot-cache"])
def test_each_committed_bundle_matches_its_manifest(name: str):
    module = _load()
    bundle = module.BUNDLES[name]
    assert bundle.manifest.is_file(), f"data/{name}.manifest.json is missing"
    assert module.verify(bundle, quiet=True), f"run scripts/derived_refs.py pack {name}"


@pytest.mark.parametrize("name", ["derived-refs", "uniprot-cache"])
def test_every_part_stays_under_githubs_warning(name: str):
    module = _load()
    manifest = json.loads(module.BUNDLES[name].manifest.read_text())
    assert manifest["parts"], f"{name} has no parts"
    for part in manifest["parts"]:
        assert part["bytes"] <= module.PART_BYTES
        assert part["bytes"] < GITHUB_WARNING, f"{part['name']} would trip GitHub's size warning"
        assert part["bytes"] < GITHUB_LIMIT


def test_the_derived_bundle_holds_what_the_module_promises():
    module = _load()
    bundle = module.BUNDLES["derived-refs"]
    manifest = json.loads(bundle.manifest.read_text())
    listed = {entry["path"] for entry in manifest["contents"]}
    required = {item.path for item in module.DERIVED if item.required}
    assert required <= listed, required - listed
    names = _names(bundle)
    for path in required:
        assert path in names, f"{path} is not in the tarball"
    # The report is scored against these; the bundle is pointless without them.
    assert "refs/reference_ranges.json" in names
    assert "refs/taxonomic_cohort.json" in names
    assert any(n.startswith("refs/age/") for n in names)


def test_the_uniprot_bundle_is_exactly_what_the_current_panels_need():
    module = _load()
    bundle = module.BUNDLES["uniprot-cache"]
    manifest = json.loads(bundle.manifest.read_text())
    listed = {entry["path"] for entry in manifest["contents"]}
    wanted = {item.path for item in module.uniprot_items()}
    assert wanted, "no UniProt entries in the panels?"
    assert listed == wanted, (
        f"the panels changed: run scripts/derived_refs.py pack uniprot-cache "
        f"(missing {sorted(wanted - listed)[:3]}, stale {sorted(listed - wanted)[:3]})"
    )
    # Every table comes with the sidecar that records its query.
    tables = {p for p in listed if p.endswith(".tsv.gz")}
    for table in tables:
        assert table[:-len(".tsv.gz")] + ".tsv.meta.json" in listed, table


def test_no_bundle_contains_sequencing_reads_or_results():
    module = _load()
    for name in module.DEFAULT_ORDER:
        for member in _names(module.BUNDLES[name]):
            assert member.startswith("refs/"), member
            assert not member.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz", ".bam", ".sam")), member
            assert ".." not in Path(member).parts, member


def test_every_listed_item_has_a_reason():
    module = _load()
    for item in module.DERIVED:
        assert item.path.startswith("refs/"), item.path
        assert len(item.why) > 20, f"{item.path} needs a real explanation"


def _names(bundle) -> set[str]:
    parts = bundle.parts()
    data = b"".join(p.read_bytes() for p in parts)
    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as tar:
        return set(tar.getnames())


# --------------------------------------------------------------------------- #
# pack / restore / verify / check
# --------------------------------------------------------------------------- #


def test_pack_then_restore_round_trips(fake_repo):
    module, root, bundle = fake_repo
    assert module.pack(bundle, quiet=True) == 0
    assert module.verify(bundle, quiet=True)
    assert module.check(bundle, quiet=True)
    manifest = json.loads(bundle.manifest.read_text())
    assert manifest["n_files"] == 4
    assert [c["path"] for c in manifest["contents"]] == [
        "refs/reference_ranges.json", "refs/age", "refs/crosswalk", "refs/panels",
    ]  # bacdive is optional and absent, so it is simply not listed

    shutil.rmtree(root / "refs")
    assert not module.check(bundle, quiet=True)
    assert module.restore(bundle, quiet=True) == 0
    assert module.check(bundle, quiet=True)
    assert (root / "refs" / "age" / "model.json").read_text() == '{"model": 1}'
    assert (root / "refs" / "panels" / "mechanisms" / "mech.bft1.fasta").read_text() == ">x\nACGT\n"
    assert (root / "refs" / "reference_ranges.json").read_text() == '{"ranges": []}'


def test_restore_leaves_existing_files_alone_unless_forced(fake_repo):
    module, root, bundle = fake_repo
    assert module.pack(bundle, quiet=True) == 0
    target = root / "refs" / "age" / "model.json"
    target.write_text("newer, local")
    assert module.restore(bundle, quiet=True) == 0
    assert target.read_text() == "newer, local"
    assert module.restore(bundle, force=True, quiet=True) == 0
    assert target.read_text() == '{"model": 1}'


def test_pack_refuses_when_a_required_item_is_missing(fake_repo, capsys):
    module, root, bundle = fake_repo
    (root / "refs" / "reference_ranges.json").unlink()
    assert module.pack(bundle, quiet=True) == 1
    assert "refs/reference_ranges.json" in capsys.readouterr().err


def test_a_tampered_archive_is_refused(fake_repo, capsys):
    module, root, bundle = fake_repo
    assert module.pack(bundle, quiet=True) == 0
    with bundle.archive.open("ab") as handle:
        handle.write(b"\0")
    assert not module.verify(bundle, quiet=True)
    assert module.restore(bundle, quiet=True) == 1
    assert "does not match" in capsys.readouterr().err


def test_packing_is_reproducible(fake_repo):
    module, root, bundle = fake_repo
    assert module.pack(bundle, quiet=True) == 0
    first = json.loads(bundle.manifest.read_text())["parts"][0]["sha256"]
    assert module.pack(bundle, quiet=True) == 0
    second = json.loads(bundle.manifest.read_text())["parts"][0]["sha256"]
    assert first == second, "the same data must pack to the same bytes"


def test_a_large_archive_is_split_and_joined(fake_repo, monkeypatch):
    module, root, bundle = fake_repo
    (root / "refs" / "age" / "big.bin").write_bytes(os.urandom(50_000))  # incompressible
    monkeypatch.setattr(module, "PART_BYTES", 16_000)
    assert module.pack(bundle, quiet=True) == 0
    parts = bundle.parts()
    assert len(parts) >= 3
    assert all(p.stat().st_size <= 16_000 for p in parts)
    assert not bundle.archive.exists()
    assert module.verify(bundle, quiet=True)
    shutil.rmtree(root / "refs")
    assert module.restore(bundle, quiet=True) == 0
    assert (root / "refs" / "age" / "big.bin").stat().st_size == 50_000


def test_restore_refuses_paths_outside_the_repository(fake_repo, tmp_path, capsys):
    module, root, bundle = fake_repo
    assert module.pack(bundle, quiet=True) == 0
    evil = io.BytesIO()
    with tarfile.open(fileobj=evil, mode="w:gz") as tar:
        info = tarfile.TarInfo("refs/../../escaped.txt")
        payload = b"nope"
        info.size = len(payload)
        tar.addfile(info, io.BytesIO(payload))
    bundle.archive.write_bytes(evil.getvalue())
    manifest = json.loads(bundle.manifest.read_text())
    manifest["parts"][0]["bytes"] = bundle.archive.stat().st_size
    manifest["parts"][0]["sha256"] = module._sha256(bundle.archive)
    bundle.manifest.write_text(json.dumps(manifest))
    assert module.restore(bundle, quiet=True) == 1
    assert "outside the repository" in capsys.readouterr().err
    assert not (tmp_path / "escaped.txt").exists()


def test_check_reports_what_is_missing(fake_repo, capsys):
    module, root, bundle = fake_repo
    assert module.pack(bundle, quiet=True) == 0
    (root / "refs" / "age").rename(root / "refs" / "age-moved")
    assert module.check(bundle) is False
    out = capsys.readouterr().out
    assert "missing from refs/" in out and "refs/age" in out


def test_cli_dispatches_every_verb(fake_repo, capsys):
    module, root, bundle = fake_repo
    assert module.main(["pack", "--quiet"]) == 0
    assert module.main(["verify", "derived-refs", "--quiet"]) == 0
    assert module.main(["check", "--quiet"]) == 0
    assert module.main(["list"]) == 0
    assert "derived-refs" in capsys.readouterr().out
    with pytest.raises(SystemExit, match="unknown bundle"):
        module.main(["verify", "nope"])
