"""The lane engines' parsers, units and floors, on fixture files.

No tools run here. Each parser is given a small file in the tool's own
format and must produce observations with the lane's own unit, a native
identifier, and the status the lane itself is allowed to give.
"""

from __future__ import annotations

import json
from pathlib import Path

from openbiota.completeness import REQUIRED, audit
from openbiota.engines import kraken, lanebase, metaphlan_jan26, motus, singlem
from openbiota.refsources import SOURCES, by_source

REPO = Path(__file__).resolve().parent.parent


def test_every_reference_family_r1_to_r10_is_declared_with_a_licence() -> None:
    srcs = by_source()
    for needed in ("mpa_jan26", "gtdb_r232", "globdb_r232", "motus_db", "hrgm2", "uhgg_v2.0.2", "uhgg_kraken2",
                   "elgg", "humgut2", "hrom", "singlem_globdb", "sgb_bridges"):
        assert needed in srcs, f"{needed} is not declared in refsources"
    for f in SOURCES:
        assert f.license, f"{f.dest} has no licence recorded"
        assert f.url.startswith(("https://", "http://")) and "latest" not in f.url.lower(), f"{f.dest}: moving alias or bad URL"


def test_mpa_jan26_parser_reads_sgb_rows_and_unclassified(tmp_path: Path) -> None:
    prof = tmp_path / "p.tsv"
    prof.write_text(
        "#mpa_vJan26_CHOCOPhlAnSGB_202605\n#/x/metaphlan ...\n#12345 reads processed\n#SampleID\tMetaphlan_Analysis\n"
        "#clade_name\tNCBI_tax_id\trelative_abundance\tadditional_species\n"
        "UNCLASSIFIED\t-1\t9.34\t\n"
        "k__Bacteria\t2\t90.66\t\n"
        "k__Bacteria|p__Bacillota|c__Clostridia|o__Lachnospirales|f__Lachnospiraceae|g__Blautia|s__Blautia_wexlerae\t1796646\t7.07\t\n"
        "k__Bacteria|p__Bacillota|c__Clostridia|o__Lachnospirales|f__Lachnospiraceae|g__Blautia|s__Blautia_wexlerae|t__SGB4837\t\t7.07\t\n"
        "k__Bacteria|p__Bacillota|c__Clostridia|o__Lachnospirales|f__Lachnospiraceae|g__GGB6612|s__GGB6612_SGB9346|t__SGB9346\t\t2.58\t\n"
    )
    obs, summary = metaphlan_jan26._parse(prof)
    assert summary["unclassified_percent"] == 9.34 and summary["reads_processed"] == 12345
    assert [o.native_id for o in obs] == ["SGB4837", "SGB9346"]
    assert obs[0].species == "Blautia wexlerae" and obs[0].abundance_unit.startswith("percent of classified")
    assert obs[1].support_metrics["unnamed_sgb"] is True
    assert summary["phyla"] == {"Bacillota": 9.65}


def test_motus_parser_uses_assigned_inserts_as_the_denominator(tmp_path: Path) -> None:
    prof = tmp_path / "m.tsv"
    prof.write_text(
        "#tool_version=4.1.0\tdatabase=4.1\nmOTU\tTaxonomy\tSAMPLE\n"
        "mOTUv4.0_000158\td__Bacteria;p__Bacillota;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Blautia_A;s__Blautia_A wexlerae\t300.0\n"
        "mOTUv4.0_000113\td__Bacteria;p__Bacillota;c__Clostridia;o__Lachnospirales;f__Lachnospiraceae;g__Lachnospira;s__Lachnospira eligens_A\t100.0\n"
        "unassigned\tunassigned\t100.0\n"
    )
    obs, summary = motus._parse(prof)
    assert len(obs) == 2 and summary["unassigned_fraction"] == 0.2
    assert obs[0].abundance_value == 75.0 and obs[0].abundance_unit == "percent of mOTU-assigned inserts"
    assert obs[0].native_id == "mOTUv4.0_000158" and obs[0].species == "Blautia A wexlerae"


def test_kraken_inventory_floor_and_supported_floor() -> None:
    names = {1: "root", 2: "Bacteria", 10: "Blautia", 11: "Blautia wexlerae", 12: "Blautia obeum"}
    nodes = {1: (1, "root"), 2: (1, "domain"), 10: (2, "genus"), 11: (10, "species"), 12: (10, "species")}
    reps = {11: "MGYG000000001", 12: "MGYG000000002"}
    report = Path("/tmp/_ob_kraken_report.tsv")
    report.write_text(
        "10.00\t100000\t100000\tU\t0\tunclassified\n"
        "90.00\t900000\t0\tR\t1\troot\n"
        "90.00\t900000\t0\tD\t2\tBacteria\n"
        "90.00\t900000\t0\tG\t10\tBlautia\n"
        "80.00\t800000\t800000\tS\t11\tBlautia wexlerae\n"
        "0.01\t120\t120\tS\t12\tBlautia obeum\n"
    )
    minimizers = {11: (250000, 900000), 12: (900, 1200)}
    obs, summary = kraken._parse(report, None, minimizers, names, nodes, reps)
    by = {o.species: o for o in obs}
    assert by["Blautia wexlerae"].status == "supported"
    assert by["Blautia wexlerae"].support_metrics["below_inventory_floor"] is False
    assert by["Blautia obeum"].status == "provisional"
    assert by["Blautia obeum"].support_metrics["below_inventory_floor"] is True, "120 reads and 900 minimizers is spillover"
    assert summary["n_species_above_inventory_floor"] == 1
    assert by["Blautia wexlerae"].support_metrics["uhgg_species_rep"] == "MGYG000000001"


def test_singlem_parser_marks_unresolved_lineages_and_triggers(tmp_path: Path) -> None:
    prof = tmp_path / "s.tsv"
    prof.write_text(
        "sample\tcoverage\ttaxonomy\n"
        "S\t24.6\tRoot; d__Bacteria; p__Bacillota; c__Clostridia; o__Lachnospirales; f__Lachnospiraceae; g__Blautia_A; s__Blautia_A wexlerae\n"
        "S\t13.4\tRoot; d__Bacteria; p__Bacillota; c__Clostridia; o__Lachnospirales; f__Lachnospiraceae\n"
        "S\t0.4\tRoot; d__Bacteria; p__Bacillota\n"
    )
    obs, summary = singlem._parse(prof)
    assert summary["n_species_resolved"] == 1
    fam = next(o for o in obs if o.rank == "family")
    assert fam.status == "ambiguous" and fam.support_metrics["assembly_trigger"] is True
    assert summary["assembly_triggers"][0]["rank"] == "family" and summary["assembly_triggers"][0]["coverage"] == 13.4
    assert all(o.abundance_unit.startswith("SingleM mean marker coverage") for o in obs)


def test_lane_result_round_trip() -> None:
    r = lanebase.LaneResult(lane_id="x", tool="t", tool_version="1", reference_release_id="rel", taxonomy_release="tax",
                            observations=[lanebase.Observation("id1", "d__B;g__G;s__G sp", "species", "supported", 1.5, "percent", "reads",
                                                                {"k": 1}, ("alt",))],
                            elapsed_s=1.0, cached=False, command=("a", "b"), input_sha256="abc", raw_result_uri="/x")
    back = lanebase.LaneResult.from_json(json.loads(json.dumps(r.to_json())))
    assert back.observations[0].alternatives == ("alt",) and back.species_percent == {"G sp": 1.5}


def test_the_gate_requires_every_new_lane_and_stage() -> None:
    keys = {s.key for s in REQUIRED}
    for k in ("lane_jan26", "lane_globdb", "lane_motus", "lane_kraken", "lane_singlem", "confirmation", "strain_analysis", "bacdive", "assembly"):
        assert k in keys
    # A lane that did not run is MISSING; a lane that ran is fine; not_triggered assembly is complete.
    r = {"detection": {"lane_status": {"metaphlan_jan26": {"state": "ran", "n_species": 250},
                                       "sylph_globdb": {"state": "not_assessed", "reason": "db missing"}}},
         "assembly": {"status": "not_triggered"}}
    rep = audit(r)
    by = {s["key"]: s for s in rep.stages}
    assert by["lane_jan26"]["ran"] and not by["lane_globdb"]["ran"] and by["assembly"]["ran"]


def test_makefile_installs_and_locks_everything() -> None:
    mk = (REPO / "Makefile").read_text()
    for target in ("tools-expanded:", "refs-expanded:", "expansion:", "expansion-crosswalks:", "expansion-reconcile:"):
        assert target in mk
    for tool in ("4.2.6", "mOTUs/archive/refs/tags/4.1.0", "singlem==0.21.4", "kraken2 bracken skani", "SPAdes", "lock_tools.py"):
        assert tool in mk, f"Makefile does not install/lock {tool}"


def test_external_html_links_open_safely() -> None:
    """Spec 0.8.4 §7: external HTML links carry target=_blank rel=noopener noreferrer.

    The report is a PDF; the only HTML the project writes today has no
    external anchors. This guards whatever HTML is written next.
    """
    import re

    offenders = []
    for path in (REPO / "openbiota").rglob("*.py"):
        for i, line in enumerate(path.read_text().splitlines(), 1):
            if re.search(r"<a\s+[^>]*href=[\"']https?://", line) and not (
                'target="_blank"' in line and 'rel="noopener noreferrer"' in line
            ):
                offenders.append(f"{path.name}:{i}")
    assert not offenders, offenders
