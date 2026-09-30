"""Profile schema and loader.

The loader's job is to refuse profiles that cannot be audited, so most of
these tests assert that something is rejected.
"""

from __future__ import annotations

import pytest
import yaml

from openbiota.errors import PanelError
from openbiota.profiles import (
    CONFLICTED_TAXA,
    DEFAULT_MODULE_WEIGHTS,
    load_profile_set,
    parse_profile,
)

MINIMAL = {
    "name": "testprofile",
    "version": "0.1.0",
    "label": "Test pattern",
    "status": "research_only",
    "summary": "A synthetic profile for tests.",
    "citation": "Nobody et al. 2026",
    "references": {
        "primary": {"source": None, "note": "none exists"},
        "secondary": {"source": "curatedMetagenomicData", "match_on": ["country", "sex"]},
    },
    "modules": {
        "taxonomic": {
            "weight_in_combined": 0.5,
            "features": [
                {
                    "name": "Faecalibacterium_prausnitzii",
                    "level": "species",
                    "direction": "decreased",
                    "w": 1.0, "q": 0.8, "s": 0.5, "c": 1.0,
                    "sources": ["Somebody 2023"],
                }
            ],
        },
        "functional": {
            "weight_in_combined": 0.35,
            "features": [
                {"name": "butyrate", "engine": "diamond", "direction": "decreased",
                 "w": 1.0, "q": 1.0, "s": 0.5, "c": 1.0}
            ],
        },
        "ecological": {
            "weight_in_combined": 0.15,
            "features": [
                {"name": "shannon_diversity", "direction": "decreased",
                 "w": 0.5, "q": 0.6, "s": 0.2, "c": 0.6}
            ],
        },
    },
}


def _profile(**overrides):
    import copy

    data = copy.deepcopy(MINIMAL)
    data.update(overrides)
    return parse_profile(data)


def test_minimal_profile_parses() -> None:
    profile = _profile()
    assert profile.name == "testprofile"
    assert len(profile.features()) == 3
    assert profile.combined_weights == {"taxonomic": 0.5, "functional": 0.35, "ecological": 0.15}


def test_direction_sign() -> None:
    profile = _profile()
    decreased = profile.features(module="taxonomic")[0]
    assert decreased.d == -1
    increased = parse_profile(
        {
            **MINIMAL,
            "modules": {
                **MINIMAL["modules"],
                "taxonomic": {
                    "weight_in_combined": 0.5,
                    "features": [
                        {"name": "Fusobacterium_nucleatum", "level": "species",
                         "direction": "increased", "w": 1.0, "q": 1.0, "s": 1.0, "c": 1.0}
                    ],
                },
            },
        }
    )
    assert increased.features(module="taxonomic")[0].d == 1


def test_scoring_weight_is_w_alone_and_audit_factor_is_the_product() -> None:
    """Evidence quality never multiplies the score (spec 7.1).

    Multiplying q, s and c into the weight pulls a weak-evidence feature
    toward the midpoint, where it reads as biologically neutral rather than
    unknown. The product is preserved for audit only.
    """
    feature = _profile().features(module="taxonomic")[0]
    assert feature.factor == pytest.approx(1.0)
    assert feature.audit_factor == pytest.approx(1.0 * 0.8 * 0.5 * 1.0)


def test_module_weights_must_sum_to_one() -> None:
    broken = {
        **MINIMAL,
        "modules": {
            "taxonomic": {"weight_in_combined": 0.9, **{"features": MINIMAL["modules"]["taxonomic"]["features"]}},
            "functional": {"weight_in_combined": 0.9, "features": MINIMAL["modules"]["functional"]["features"]},
        },
    }
    with pytest.raises(PanelError, match="must sum to 1.0"):
        parse_profile(broken)


def test_version_is_required_and_semantic() -> None:
    with pytest.raises(PanelError, match="semantic version"):
        _profile(version="draft")


def test_taxonomic_feature_must_declare_level() -> None:
    data = {
        **MINIMAL,
        "modules": {
            **MINIMAL["modules"],
            "taxonomic": {
                "weight_in_combined": 0.5,
                "features": [{"name": "Some_species", "direction": "decreased"}],
            },
        },
    }
    with pytest.raises(PanelError, match="must declare 'level'"):
        parse_profile(data)


def _with_taxon(name: str, level: str):
    return {
        **MINIMAL,
        "modules": {
            **MINIMAL["modules"],
            "taxonomic": {
                "weight_in_combined": 0.5,
                "features": [{"name": name, "level": level, "direction": "decreased"}],
            },
        },
    }


@pytest.mark.parametrize("taxon", sorted(CONFLICTED_TAXA))
def test_conflicted_taxon_is_refused_at_any_level(taxon: str) -> None:
    """Aggregating a conflicted genus reverses direction in practice.

    Refused whichever level is declared, including a bare genus name
    mislabelled as species.
    """
    for level in ("genus", "species"):
        with pytest.raises(PanelError, match="conflicted-taxa list"):
            parse_profile(_with_taxon(taxon.capitalize(), level))


def test_genus_level_is_never_accepted_for_shotgun_evidence() -> None:
    with pytest.raises(PanelError, match="Genus level is never accepted"):
        parse_profile(_with_taxon("Akkermansia", "genus"))


def test_genus_level_binds_only_by_declared_amplicon_transport() -> None:
    """A 16S genus finding may bind, but only labelled as a transport (spec 5.2)."""
    data = {
        **MINIMAL,
        "spec_version": 3,
        "modules": {
            "lee_2024": {
                "evidence_type": "TAX_REL",
                "assay_transport": "amplicon_taxon_to_shotgun",
                "features": [{"name": "Akkermansia", "level": "genus", "direction": "increased"}],
            },
        },
    }
    profile = parse_profile(data)
    module = profile.modules["lee_2024"]
    assert module.assay_transport == "amplicon_taxon_to_shotgun"
    assert module.features[0].level == "genus"
    assert profile.is_transported


def test_conflicted_genus_binds_only_with_material_conflict_declared() -> None:
    data = {
        **MINIMAL,
        "spec_version": 3,
        "modules": {
            "lee_2024": {
                "evidence_type": "TAX_REL",
                "assay_transport": "amplicon_taxon_to_shotgun",
                "features": [
                    {"name": "Bacteroides", "level": "genus", "direction": "decreased",
                     "direction_conflict": "material"},
                ],
            },
        },
    }
    profile = parse_profile(data)
    assert profile.has_direction_conflict
    data["modules"]["lee_2024"]["features"][0].pop("direction_conflict")
    with pytest.raises(PanelError, match="conflicted-taxa list"):
        parse_profile(data)


def test_unbound_evidence_type_loads_and_says_what_would_bind_it() -> None:
    """A KO panel loads without an engine and reports ``would_bind_if`` (spec 13)."""
    data = {
        **MINIMAL,
        "spec_version": 3,
        "modules": {
            **{k: {**v, "evidence_type": {"taxonomic": "TAX_REL", "functional": "GENE_ABUND",
                                          "ecological": "ECO", "phenotype": "ECO"}[k]}
               for k, v in MINIMAL["modules"].items()},
            "metaad_ko": {
                "evidence_type": "PATH_ABUND",
                "features": [{"name": "K00857", "direction": "increased"}],
            },
        },
    }
    profile = parse_profile(data)
    ko = profile.modules["metaad_ko"]
    assert not ko.bound
    assert not ko.fuses
    assert "HUMAnN" in (ko.would_bind_if or "")
    assert "metaad_ko" not in profile.combined_weights


def test_same_study_group_shares_one_vote() -> None:
    """Two modules from one cohort's stool DNA are one evidence group (spec 9.1)."""
    data = {
        **MINIMAL,
        "spec_version": 3,
        "modules": {
            "taxa_a": {"evidence_type": "TAX_REL", "study_group": "cohort_a",
                       "features": [{"name": "Species_one", "level": "species", "direction": "increased"}]},
            "genes_a": {"evidence_type": "GENE_ABUND", "study_group": "cohort_a",
                        "features": [{"name": "butyrate", "direction": "decreased"}]},
            "taxa_b": {"evidence_type": "TAX_REL",
                       "features": [{"name": "Species_two", "level": "species", "direction": "increased"}]},
        },
    }
    profile = parse_profile(data)
    weights = profile.combined_weights
    assert weights["taxa_a"] == pytest.approx(0.25)
    assert weights["genes_a"] == pytest.approx(0.25)
    assert weights["taxa_b"] == pytest.approx(0.5)


def test_sensitivity_panel_never_fuses() -> None:
    data = {
        **MINIMAL,
        "spec_version": 3,
        "modules": {
            "strict": {"evidence_type": "TAX_REL",
                       "features": [{"name": "Species_one", "level": "species", "direction": "increased"}]},
            "relaxed": {"evidence_type": "TAX_REL", "type": "sensitivity_panel",
                        "features": [{"name": "Species_two", "level": "species", "direction": "increased"}]},
        },
    }
    profile = parse_profile(data)
    assert profile.modules["relaxed"].status == "exploratory"
    assert list(profile.combined_weights) == ["strict"]


def test_carrier_abundance_cannot_be_declared_observed() -> None:
    data = {
        **MINIMAL,
        "spec_version": 3,
        "modules": {
            "carriers": {"evidence_type": "CARRIER_ABUNDANCE", "claim_level": "observed",
                         "features": [{"name": "butyrate", "direction": "decreased"}]},
        },
    }
    with pytest.raises(PanelError, match="carrier_proxy"):
        parse_profile(data)


def test_conflicted_species_is_allowed() -> None:
    """The species is fine; only the genus is refused."""
    profile = _profile()
    assert profile.features(module="taxonomic")[0].name == "Faecalibacterium_prausnitzii"


def test_weight_factors_must_be_unit_interval() -> None:
    data = {
        **MINIMAL,
        "modules": {
            **MINIMAL["modules"],
            "taxonomic": {
                "weight_in_combined": 0.5,
                "features": [
                    {"name": "X_y", "level": "species", "direction": "decreased", "q": 1.4}
                ],
            },
        },
    }
    with pytest.raises(PanelError, match=r"must be a number in \[0, 1\]"):
        parse_profile(data)


def test_phenotype_module_is_excluded_from_combined() -> None:
    data = {
        **MINIMAL,
        "modules": {
            **MINIMAL["modules"],
            "phenotype": {
                "weight_in_combined": 0.4,  # should be forced to zero
                "features": [
                    {"name": "A_b", "level": "species", "direction": "increased", "w": 0.25}
                ],
            },
        },
    }
    profile = parse_profile(data)
    assert profile.modules["phenotype"].weight_in_combined == 0.0
    assert "phenotype" not in profile.combined_weights


def test_duration_dependent_requires_strata() -> None:
    with pytest.raises(PanelError, match="at least two strata"):
        _profile(duration_dependent=True)


def test_cross_engine_check_must_reference_real_features() -> None:
    data = {
        **MINIMAL,
        "cross_engine_checks": [
            {"taxonomic_feature": "Nonexistent_species", "functional_feature": "butyrate"}
        ],
    }
    with pytest.raises(PanelError, match="is not a feature of this profile"):
        parse_profile(data)


def test_abstention_rules_get_default_reasons() -> None:
    profile = _profile(abstain_if=[{"usable_nonhost_reads_lt": 500000}])
    assert profile.abstain_if[0].kind == "usable_nonhost_reads_lt"
    assert "500,000" in profile.abstain_if[0].reason


def test_default_module_weights_are_the_documented_split() -> None:
    assert DEFAULT_MODULE_WEIGHTS == {
        "taxonomic": 0.50, "functional": 0.35, "ecological": 0.15
    }


# --------------------------------------------------------------------------- #
# the shipped profiles
# --------------------------------------------------------------------------- #


@pytest.fixture(scope="module")
def shipped():
    from pathlib import Path

    return load_profile_set(Path(__file__).resolve().parents[1] / "profiles")


def test_shipped_profiles_load(shipped) -> None:
    names = {p.name for p in shipped.profiles}
    # The original three, plus the spec-3 library.
    assert {"crc", "mecfs", "longcovid"} <= names
    assert {
        "ad_clinical", "ad_mci", "ad_preclinical_amyloid", "alopecia_areata",
        "androgenetic_alopecia", "ibs", "ibs_c", "ibs_d", "ibs_m", "crohns", "uc", "ibd",
        "t2d", "cirrhosis", "parkinsons", "ms", "acvd", "ra", "masld", "mdd", "t1d",
        "celiac", "ankylosing_spondylitis", "hypertension", "ckd", "obesity", "adenoma",
    } <= names
    assert len(names) >= 30


def test_crc_is_marked_a_validation_control(shipped) -> None:
    crc = shipped.by_name("crc")
    assert crc.status == "validation_control"
    assert any("NOT A CANCER SCREEN" in c.upper() for c in crc.caveats)


def test_mecfs_is_duration_dependent_with_two_strata(shipped) -> None:
    mecfs = shipped.by_name("mecfs")
    assert mecfs.duration_dependent
    assert {s.name for s in mecfs.strata} == {"short_term", "long_term"}


def test_mecfs_carries_the_cross_engine_check(shipped) -> None:
    mecfs = shipped.by_name("mecfs")
    assert len(mecfs.cross_engine_checks) == 1
    check = mecfs.cross_engine_checks[0]
    assert check.taxonomic_feature == "Faecalibacterium_prausnitzii"
    assert check.functional_feature == "butyrate"


def test_longcovid_penalises_the_single_institution_family(shipped) -> None:
    """Tier A features from one recruitment family get roughly one vote."""
    longcovid = shipped.by_name("longcovid")
    tier_a = [f for f in longcovid.features(module="taxonomic") if f.tier == "A"]
    assert tier_a, "expected tier A features"
    assert all(f.q <= 0.5 for f in tier_a), "q should reflect one independent vote"


def test_longcovid_declares_its_discrimination_unmeasured(shipped) -> None:
    note = shipped.by_name("longcovid").effect_size_note
    assert "UNMEASURED" in note
    assert "no AUC is printed" in note


def test_no_profile_uses_a_conflicted_genus(shipped) -> None:
    """A conflicted genus is refused for shotgun evidence outright, and for
    transported 16S evidence is admitted only when it says so (genus level,
    material direction conflict, amplicon transport)."""
    for profile in shipped.profiles:
        for module in profile.modules.values():
            for feature in module.features:
                if feature.name.strip().lower() not in CONFLICTED_TAXA:
                    continue
                assert feature.level == "genus", (profile.name, module.name, feature.name)
                assert feature.direction_conflict == "material", (profile.name, feature.name)
                assert module.assay_transport == "amplicon_taxon_to_shotgun", (
                    profile.name, module.name, feature.name,
                )


def test_every_feature_has_all_four_factors(shipped) -> None:
    for profile in shipped.profiles:
        for feature in profile.features():
            for name in ("w", "q", "s", "c"):
                assert 0.0 <= getattr(feature, name) <= 1.0


def test_shipped_yaml_is_valid_yaml() -> None:
    from pathlib import Path

    for path in (Path(__file__).resolve().parents[1] / "profiles").glob("*.yaml"):
        assert yaml.safe_load(path.read_text(encoding="utf-8"))
