"""Lane cohorts: a population per lane, one population per organism, named.

Pure fixtures. What is pinned:
  * an organism the scoring population ranked keeps that rank and is
    labelled "scoring cohort"; a lane cohort never overrides it
  * an organism seen only by a lane with its own cohort gets a percentile
    from that cohort and carries the population's name
  * an organism seen only by a lane without a cohort stays unranked
  * with lane cohorts switched off (empty mapping) nothing changes
  * the inventory reports how many organisms each population ranked
  * the JSON round trip keeps the population name
"""

from __future__ import annotations

import random

from openbiota import inventory
from openbiota.expansion import cohorts
from openbiota.refcohort import CohortSample, ReferenceCohort


def _lane_cohort(taxa: list[str], n: int = 60, seed: int = 7) -> cohorts.LaneCohort:
    rng = random.Random(seed)
    abundance = [[round(max(0.0, rng.gauss(mu, mu / 2)), 4) if rng.random() < prev else 0.0 for _ in range(n)]
                 for mu, prev in ((5.0, 0.95), (1.0, 0.8), (0.2, 0.5), (3.0, 0.9))][: len(taxa)]
    ids = [f"S{i}" for i in range(n)]
    meta = {s: CohortSample(sample_id=s, study_name="PRJTEST", condition="adult stool", country=None, age=None, sex=None,
                            bmi=None, n_reads=3_000_000, platform=None, antibiotics=None, westernized=None) for s in ids}
    ref = ReferenceCohort(snapshot="GlobDB r232", profiler="sylph 1.0.0", taxa=taxa, sample_ids=ids,
                          abundance=abundance, metadata=meta, manifest_id="t")
    return cohorts.LaneCohort(lane="globdb", cohort=ref, provenance={
        "study": "PRJTEST", "lane_label": "sylph 1.0.0 against GlobDB r232", "depth_read_pairs": 3_000_000,
        "sampling": "prefix"})


class _Ranked:
    def __init__(self, species: str, percentile: float, prevalence: float) -> None:
        self.species, self.percentile, self.cohort_prevalence = species, percentile, prevalence
        self.groups, self.trace, self.rare_in_cohort = (), False, False


def _lanes():
    scoring = inventory.Lane(name="scoring", species={"Blautia_wexlerae": 7.0}, rankable=True, method="marker")
    jan26 = inventory.Lane(name="jan26", species={"Blautia_wexlerae": 7.2, "Faecalibacterium_prausnitzii": 3.0},
                           rankable=False, primary=True, method="marker")
    globdb = inventory.Lane(name="globdb", species={"Blautia_wexlerae": 6.9, "Blautia_A_MGYG000001338": 6.35,
                                                    "Gemmiger_MOTU40_023295": 0.03}, rankable=False, method="genome_sketch")
    motus = inventory.Lane(name="motus", species={"Ruminococcus_bromii": 2.0}, rankable=False, method="universal_marker")
    return [scoring, jan26, globdb, motus]


TAXA = ["Blautia_wexlerae", "Blautia_A_MGYG000001338", "Gemmiger_MOTU40_023295", "Faecalibacterium_prausnitzii"]


def test_scoring_rank_is_kept_and_named() -> None:
    inv = inventory.build(_lanes(), ranked_rows=[_Ranked("Blautia_wexlerae", 81.0, 0.97)],
                          lane_cohorts={"globdb": _lane_cohort(TAXA)})
    o = inv.get("Blautia_wexlerae")
    assert o is not None and o.percentile == 81.0 and o.percentile_source == "scoring cohort"


def test_lane_only_organism_is_ranked_in_its_own_population() -> None:
    lc = _lane_cohort(TAXA)
    inv = inventory.build(_lanes(), ranked_rows=[_Ranked("Blautia_wexlerae", 81.0, 0.97)], lane_cohorts={"globdb": lc})
    o = inv.get("Blautia_A_MGYG000001338")
    assert o is not None
    assert o.percentile is not None and 0.0 <= o.percentile <= 100.0
    assert o.prevalence is not None and 0.0 < o.prevalence <= 1.0
    assert o.percentile_source == "globdb cohort, n=60"
    # 6.35% in a population whose carriers centre on 1% is high
    assert o.percentile > 75.0
    assert inv.percentile_sources == {"scoring cohort": 1, "globdb cohort, n=60": 2}


def test_lane_without_a_cohort_stays_unranked() -> None:
    inv = inventory.build(_lanes(), ranked_rows=[_Ranked("Blautia_wexlerae", 81.0, 0.97)],
                          lane_cohorts={"globdb": _lane_cohort(TAXA)})
    bromii = inv.get("Ruminococcus_bromii")
    assert bromii is not None and bromii.percentile is None and bromii.percentile_source is None
    # seen by the primary lane only, whose namespace has no population here
    fp = inv.get("Faecalibacterium_prausnitzii")
    assert fp is not None and fp.percentile is None


def test_switched_off_changes_nothing() -> None:
    inv = inventory.build(_lanes(), ranked_rows=[_Ranked("Blautia_wexlerae", 81.0, 0.97)], lane_cohorts={})
    assert [o.species for o in inv.rankable] == ["Blautia_wexlerae"]
    assert inv.percentile_sources == {"scoring cohort": 1}


def test_round_trip_keeps_the_population_name() -> None:
    inv = inventory.build(_lanes(), ranked_rows=[_Ranked("Blautia_wexlerae", 81.0, 0.97)],
                          lane_cohorts={"globdb": _lane_cohort(TAXA)})
    back = inventory.from_json(inv.to_json())
    assert back is not None
    assert back.get("Blautia_A_MGYG000001338").percentile_source == "globdb cohort, n=60"
    assert back.percentile_sources == inv.percentile_sources


def test_cohort_file_round_trip(tmp_path) -> None:
    """A cohort written the way the builder writes it loads through load_all."""
    import json

    lc = _lane_cohort(TAXA)
    payload = {"manifest": {"manifest_id": "t", "snapshot": "GlobDB r232", "profiler": "sylph 1.0.0", "built_at": "",
                            "n_samples": 60, "n_taxa": 4, "lane": "globdb", "study": "PRJTEST",
                            "lane_label": "sylph 1.0.0 against GlobDB r232", "depth_read_pairs": 3000000, "sampling": "prefix"},
               "taxa": lc.cohort.taxa, "sample_ids": lc.cohort.sample_ids, "abundance": lc.cohort.abundance,
               "metadata": {k: v.to_json() for k, v in lc.cohort.metadata.items()}}
    (tmp_path / "globdb.cohort.json").write_text(json.dumps(payload))
    loaded = cohorts.load_all(tmp_path)
    assert set(loaded) == {"globdb"} and loaded["globdb"].n == 60
    assert "PRJTEST" in loaded["globdb"].description and "3,000,000 read pairs" in loaded["globdb"].description
