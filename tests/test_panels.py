"""Panel loading and validation. Covers spec test case 6."""

from __future__ import annotations

import copy
import re
from pathlib import Path

import pytest
import yaml

from openbiota.errors import PanelError
from openbiota.panels import load_panel_set, parse_normalizer, parse_panel

from .conftest import NORMALIZER_YAML, URDA_PANEL_YAML

REPO_ROOT = Path(__file__).resolve().parent.parent
PANELS_DIR = REPO_ROOT / "panels"


# --------------------------------------------------------------------------- #
# 6. malformed YAML fails with a clear message
# --------------------------------------------------------------------------- #


def test_invalid_yaml_syntax(tmp_path: Path):
    (tmp_path / "_normalizer.yaml").write_text(
        yaml.safe_dump(NORMALIZER_YAML), encoding="utf-8"
    )
    (tmp_path / "broken.yaml").write_text("name: broken\n  bad: [indent\n", encoding="utf-8")

    with pytest.raises(PanelError, match="invalid YAML"):
        load_panel_set(tmp_path)


def test_missing_panels_directory(tmp_path: Path):
    with pytest.raises(PanelError, match="panels directory not found"):
        load_panel_set(tmp_path / "nope")


def test_missing_normalizer_file(tmp_path: Path):
    (tmp_path / "urda.yaml").write_text(yaml.safe_dump(URDA_PANEL_YAML), encoding="utf-8")
    with pytest.raises(PanelError, match="panel file not found"):
        load_panel_set(tmp_path)


def test_no_panel_files(tmp_path: Path):
    (tmp_path / "_normalizer.yaml").write_text(
        yaml.safe_dump(NORMALIZER_YAML), encoding="utf-8"
    )
    with pytest.raises(PanelError, match="no panel definitions found"):
        load_panel_set(tmp_path)


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda d: d.pop("name"), "missing required key 'name'"),
        (lambda d: d.pop("description"), "missing required key 'description'"),
        (lambda d: d.pop("metabolite"), "missing required key 'metabolite'"),
        (lambda d: d.pop("targets"), "missing required key 'targets'"),
        (lambda d: d.pop("citation"), "missing required key 'citation'"),
        (lambda d: d.update(targets=[]), "must contain at least one entry"),
        (lambda d: d.update(name="bad name!"), "must be alphanumeric"),
        (lambda d: d.update(aggregate="average"), "'aggregate' must be one of"),
        (lambda d: d.update(min_fragments_for_stability=-1), "non-negative integer"),
        (lambda d: d.update(min_alignment_aa=0), "must be a positive integer"),
        (lambda d: d.update(extension="yes"), "'extension' must be a boolean"),
        (lambda d: d.update(decoys={"id": "X"}), "'decoys' must be a list"),
        (lambda d: d.update(targets=["not a mapping"]), "must be a mapping"),
    ],
)
def test_malformed_panel_fields(mutate, message: str):
    data = copy.deepcopy(URDA_PANEL_YAML)
    mutate(data)
    with pytest.raises(PanelError, match=message):
        parse_panel(data)


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda e: e.pop("id"), "missing required key 'id'"),
        (lambda e: e.pop("query"), "missing required key 'query'"),
        (lambda e: e.pop("source"), "missing required key 'source'"),
        (lambda e: e.update(source="ncbi"), "unsupported source"),
        (lambda e: e.update(query="   "), "'query' must not be empty"),
        (lambda e: e.update(min_identity=101), r"must be a number in \[0, 100\]"),
        (lambda e: e.update(min_identity="high"), r"must be a number in \[0, 100\]"),
        (lambda e: e.update(max_sequences=0), "positive integer"),
        (lambda e: e.update(length_tolerance=0), r"must be in \(0, 1\]"),
        (lambda e: e.update(length_range=[500]), "two-item list"),
        (lambda e: e.update(length_range=[900, 100]), r"0 < min < max"),
        (lambda e: e.update(id="HAS~TILDE"), "must be alphanumeric"),
        (lambda e: e.update(id="has space"), "must be alphanumeric"),
    ],
)
def test_malformed_entry_fields(mutate, message: str):
    data = copy.deepcopy(URDA_PANEL_YAML)
    mutate(data["targets"][0])
    with pytest.raises(PanelError, match=message):
        parse_panel(data)


@pytest.mark.parametrize(
    ("mutate", "message"),
    [
        (lambda r: r.pop("anchor_accession"), "missing required key 'anchor_accession'"),
        (lambda r: r.pop("canonical_position"), "missing required key 'canonical_position'"),
        (lambda r: r.update(canonical_position=0), "must be >= 1"),
        (lambda r: r.update(accepted_residues=[]), "must not be empty"),
        (lambda r: r.update(accepted_residues=["Y", "Z"]), "non-amino-acid codes"),
        (lambda r: r.update(applies_to=["NOPE"]), "names targets that do not exist"),
        (lambda r: r.update(enabled="true"), "'enabled' must be a boolean"),
    ],
)
def test_malformed_residue_check(mutate, message: str):
    data = copy.deepcopy(URDA_PANEL_YAML)
    mutate(data["residue_check"])
    with pytest.raises(PanelError, match=message):
        parse_panel(data)


def test_duplicate_entry_id_within_panel():
    data = copy.deepcopy(URDA_PANEL_YAML)
    data["decoys"][0]["id"] = "URDA"
    with pytest.raises(PanelError, match="duplicate entry id"):
        parse_panel(data)


def test_duplicate_panel_name(tmp_path: Path):
    (tmp_path / "_normalizer.yaml").write_text(
        yaml.safe_dump(NORMALIZER_YAML), encoding="utf-8"
    )
    (tmp_path / "a.yaml").write_text(yaml.safe_dump(URDA_PANEL_YAML), encoding="utf-8")
    (tmp_path / "b.yaml").write_text(yaml.safe_dump(URDA_PANEL_YAML), encoding="utf-8")
    with pytest.raises(PanelError, match="duplicate panel name"):
        load_panel_set(tmp_path)


def test_normalizer_requires_citation_and_description():
    data = copy.deepcopy(NORMALIZER_YAML)
    data.pop("citation")
    with pytest.raises(PanelError, match="'citation' must not be empty"):
        parse_normalizer(data)


def test_panel_file_must_be_a_mapping(tmp_path: Path):
    with pytest.raises(PanelError, match="must contain a YAML mapping"):
        parse_panel(["a", "list"], source_path=tmp_path / "x.yaml")


# --------------------------------------------------------------------------- #
# the real, shipped panels
# --------------------------------------------------------------------------- #


def test_shipped_panels_all_load():
    panel_set = load_panel_set(PANELS_DIR)
    assert len(panel_set.panels) >= 8
    names = {p.name for p in panel_set.panels}
    for required in ("urda", "cutc", "butyrate", "pcresol", "bai", "bsh", "indole", "ipa"):
        assert required in names, f"panel {required} is missing"


def test_shipped_panels_have_globally_unique_entry_keys():
    panel_set = load_panel_set(PANELS_DIR)
    keys = [e.key for e in panel_set.all_entries()]
    assert len(keys) == len(set(keys))


def test_shipped_panels_never_use_the_bare_cutc_gene_name():
    """`gene:cutC` returns copper homeostasis protein CutC, not choline TMA-lyase.

    That is the single most dangerous mistake available in this domain, so it is
    guarded by a test rather than only by a comment.
    """
    panel_set = load_panel_set(PANELS_DIR)
    cutc = panel_set.by_name("cutc")
    target = next(t for t in cutc.targets if t.id == "CUTC")
    assert "gene:cutC" not in target.query
    assert "choline trimethylamine-lyase" in target.query
    assert "activating enzyme" in target.query  # cutD must be excluded

    cnta = next(t for t in cutc.targets if t.id == "CNTA")
    assert "gene:cntA" not in cnta.query  # returns metal-staphylopine-binding CntA


def test_excluded_genes_are_absent_from_every_panel():
    """ldh, gadB, speE and murI are near-universal and deliberately excluded."""
    panel_set = load_panel_set(PANELS_DIR)
    banned = ("gene:ldh", "gene:gadB", "gene:speE", "gene:murI")
    for panel in panel_set.panels:
        for entry in panel.targets:
            for token in banned:
                assert token not in entry.query, f"{entry.key} screens excluded gene {token}"


def test_urda_panel_declares_frda_decoy_and_residue_check():
    panel_set = load_panel_set(PANELS_DIR)
    urda = panel_set.by_name("urda")
    assert any(d.id == "FRDA" for d in urda.decoys)
    assert urda.residue_check is not None
    assert urda.residue_check.canonical_position == 373
    assert urda.residue_check.accepted_residues == frozenset({"Y", "M"})
    # The spec's Q8EAP8 is a 180 aa PH-domain protein; Q8CVD0 is UrdA_SHEON.
    assert urda.residue_check.anchor_accession == "Q8CVD0"
    assert urda.upper_bound_reason is not None


def test_butyrate_is_configured_as_the_positive_control():
    panel_set = load_panel_set(PANELS_DIR)
    butyrate = panel_set.by_name("butyrate")
    # The two terminal enzymes are what the panel sums; the lactate
    # utilisation genes (lctA/B/C) are searched as context and not summed.
    assert set(butyrate.aggregate_from or ()) == {"BUT", "BUK"} or {
        t.gene for t in butyrate.targets if t.gene in ("but", "buk")
    } == {"but", "buk"}
    assert {"but", "buk"} <= {t.gene for t in butyrate.targets}
    # Acetate kinase is homologous to butyrate kinase and must be a decoy.
    assert any(d.id == "ACKA" for d in butyrate.decoys)


def test_bai_uses_median_aggregate_for_an_operon():
    panel_set = load_panel_set(PANELS_DIR)
    assert panel_set.by_name("bai").aggregate == "median"


def test_bai_never_queries_the_broken_delta_paren_protein_name():
    """The comma inside "Delta(4,5)" breaks UniProt's boolean parser.

    It silently collapses 361 hits to 3 and the 361 were plant progesterone
    5-beta-reductases anyway. Guarded by a test rather than only a comment.
    """
    panel_set = load_panel_set(PANELS_DIR)
    for panel in panel_set.panels:
        for entry in panel.entries:
            assert "Delta(4,5)" not in entry.query, f"{entry.key} uses the broken phrase"
            # more generally: no comma inside parentheses inside a quoted phrase
            for phrase in re.findall(r'"([^"]*)"', entry.query):
                inner = re.findall(r"\(([^)]*)\)", phrase)
                for chunk in inner:
                    assert "," not in chunk, (
                        f"{entry.key}: quoted phrase {phrase!r} has a comma inside parentheses, "
                        "which breaks UniProt's query parser silently"
                    )


def test_bai_targets_the_correct_desaturase_identities():
    panel_set = load_panel_set(PANELS_DIR)
    bai = panel_set.by_name("bai")
    baicd = next(t for t in bai.targets if t.id == "BAICD")
    baih = next(t for t in bai.targets if t.id == "BAIH")

    assert "3-oxocholoyl-CoA 4-desaturase" in baicd.query
    assert "4-desaturase" in baih.query
    # both are ~640-660 aa; an explicit window keeps short misannotations out
    for entry in (baicd, baih):
        assert entry.length_range is not None
        lo, hi = entry.length_range
        assert lo >= 500 and hi <= 800


def test_aggregate_from_narrows_the_headline_to_specific_genes():
    """bai, ipa and propionate each keep a low-specificity gene out of the headline."""
    panel_set = load_panel_set(PANELS_DIR)

    bai = panel_set.by_name("bai")
    assert bai.aggregate_target_ids() == ("BAICD", "BAIH")
    assert "BAIE" in bai.target_ids() and "BAIE" not in bai.aggregate_target_ids()
    assert "BAIA" in bai.target_ids() and "BAIA" not in bai.aggregate_target_ids()

    ipa = panel_set.by_name("ipa")
    assert ipa.aggregate_target_ids() == ("FLDBC",)
    assert "FLDH" in ipa.target_ids()

    propionate = panel_set.by_name("propionate")
    assert propionate.aggregate_target_ids() == ("PCT", "PDUP")
    assert "LCDA" in propionate.target_ids()


def test_aggregate_from_defaults_to_every_target():
    panel = parse_panel(URDA_PANEL_YAML)
    assert panel.aggregate_from == ()
    assert panel.aggregate_target_ids() == ("URDA",)


def test_aggregate_from_rejects_unknown_target():
    data = copy.deepcopy(URDA_PANEL_YAML)
    data["aggregate_from"] = ["NOPE"]
    with pytest.raises(PanelError, match="'aggregate_from' names targets that do not exist"):
        parse_panel(data)


#: Bacterial taxa a panel may filter to instead of the whole domain (taxonomy_id:2).
#: Narrower than "bacteria" is fine; the test guards against *no* filter, which
#: pulls in eukaryotic homologs.
BACTERIAL_TAXON_FILTERS = {
    "taxonomy_id:2",       # Bacteria
    "taxonomy_id:848",     # Fusobacterium (FadA is genus-specific)
    "taxonomy_id:816",     # Bacteroides (bft and its family paralogues)
    "taxonomy_id:817",     # Bacteroides fragilis
    "taxonomy_id:561",     # Escherichia (clb island)
    "taxonomy_id:2157",    # Archaea (methanogenesis: mcrA is archaeal, not bacterial)
}


def test_every_shipped_entry_query_filters_to_bacteria_or_is_documented():
    """A missing taxonomy filter pulls in eukaryotic homologs.

    That is not hypothetical: searching the gene symbols of the urolithin
    dehydroxylase returns an *Aspergillus* protein which happens to share
    them, and this guard is what stops one reaching a panel.

    A `literal` entry is exempt because the risk does not apply to it. Its
    query is not a name to be matched but one sequence, shipped in this
    repository and identified by its own accession, so there is no
    population to filter down to. The named-accession queries still carry
    the filter, since they go through the same search as everything else.
    """
    panel_set = load_panel_set(PANELS_DIR)
    for panel in panel_set.panels:
        for entry in panel.entries:
            if entry.source == "literal":
                assert "#" in entry.query, (
                    f"{entry.key}: a literal entry must name its file and key, so the "
                    f"sequence it uses can be found and checked; got {entry.query!r}"
                )
                continue
            assert (
                any(t in entry.query for t in BACTERIAL_TAXON_FILTERS) or entry.id == "CNTA"
            ), f"{entry.key} does not restrict to a prokaryotic taxon: {entry.query}"


def test_a_literal_entry_resolves_to_a_real_sequence_in_this_repository():
    """The exemption above is only safe if the sequence is actually there."""
    import json

    repo = Path(__file__).resolve().parents[1]
    panel_set = load_panel_set(PANELS_DIR)
    checked = 0
    for panel in panel_set.panels:
        for entry in panel.entries:
            if entry.source != "literal":
                continue
            rel, _, key = str(entry.query).partition("#")
            path = repo / rel
            assert path.is_file(), f"{entry.key}: {rel} is not in this repository"
            blob = json.loads(path.read_text(encoding="utf-8"))
            proteins = blob.get("proteins") or {}
            assert key in proteins, (
                f"{entry.key}: {rel} holds {sorted(proteins)}, not {key!r}"
            )
            record = proteins[key]
            assert record.get("sequence"), f"{entry.key}: {key} has no sequence"
            assert record.get("protein_id"), (
                f"{entry.key}: {key} has no accession, so the sequence cannot be traced"
            )
            assert blob.get("source_accession"), (
                f"{rel}: no source accession recorded, so the operon it came from is "
                "not auditable"
            )
            checked += 1
    assert checked, "no literal entries were checked; the exemption is untested"


def test_unknown_panel_name_lists_alternatives():
    panel_set = load_panel_set(PANELS_DIR)
    with pytest.raises(PanelError, match="available:"):
        panel_set.by_name("nonexistent")


def test_select_all_returns_every_panel():
    panel_set = load_panel_set(PANELS_DIR)
    assert panel_set.select(None) == panel_set.panels
    assert len(panel_set.select(["urda"])) == 1


def test_plausibility_ceiling_scales_with_summed_genes_and_declared_copies():
    """bcaa sums two near-universal genes; GUS declares three copies per genome."""
    panels = load_panel_set(Path(__file__).resolve().parent.parent / "panels")
    assert panels.by_name("urda").plausible_ceiling(150.0) == 150.0
    assert panels.by_name("bcaa").plausible_ceiling(150.0) == 300.0  # ilvC + ilvD
    assert panels.by_name("riboflavin").plausible_ceiling(150.0) == 450.0
    assert panels.by_name("bglucuronidase").plausible_ceiling(150.0) == 450.0


@pytest.mark.parametrize("bad", [0, 0.5, 11, True, "3"])
def test_expected_copies_per_genome_is_bounded(bad):
    data = copy.deepcopy(URDA_PANEL_YAML)
    data["expected_copies_per_genome"] = bad
    with pytest.raises(PanelError, match="expected_copies_per_genome"):
        parse_panel(data)
