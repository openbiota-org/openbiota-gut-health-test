"""The dietary substrate panel (A01) and its specificity rules.

The specification names four overclaims by hand. Each has a test, because
each is a thing this module exists to refuse.
"""

from __future__ import annotations

import pytest

from openbiota.extension import substrates as S


def _hits(*pairs: tuple[str, str], fragments: int = 100) -> list[S.GeneHit]:
    return [S.GeneHit(gene, family, fragments) for gene, family in pairs]


# --------------------------------------------------------------------------- #
# the fourteen views
# --------------------------------------------------------------------------- #


def test_all_fourteen_substrates_are_present_and_distinct():
    assert len(S.SUBSTRATES) == 14
    ids = [s.substrate_id for s in S.SUBSTRATES]
    assert len(set(ids)) == 14
    assert all(i.startswith("carb.") for i in ids)
    expected = {
        "cellulose", "starch_general", "resistant_starch", "chitin", "pectin", "inulin",
        "fos", "gos", "xos", "imo", "lactose", "beta_glucan", "arabinoxylan",
        "galactomannan",
    }
    assert {i.split(".", 1)[1] for i in ids} == expected


def test_every_substrate_can_be_told_apart_from_its_class():
    """A view with nothing discriminating would read as specific by default."""
    for substrate in S.SUBSTRATES:
        assert substrate.discriminating_families or substrate.discriminating_context, (
            f"{substrate.substrate_id} has nothing that separates it from its class"
        )


def test_every_shared_substrate_says_what_the_sharing_does_not_prove():
    for substrate in S.SUBSTRATES:
        if substrate.shares_with:
            assert substrate.must_not_conclude, substrate.substrate_id
            for other in substrate.shares_with:
                assert other in S.BY_ID, f"{substrate.substrate_id} shares with unknown {other}"


def test_sharing_is_declared_from_both_sides():
    """If inulin shares GH32 with FOS, FOS says so too."""
    for substrate in S.SUBSTRATES:
        for other_id in substrate.shares_with:
            other = S.BY_ID[other_id]
            assert substrate.substrate_id in other.shares_with, (
                f"{substrate.substrate_id} declares sharing with {other_id}, "
                f"but {other_id} does not declare it back"
            )


def test_a_substrate_without_a_required_family_is_refused():
    with pytest.raises(S.SubstrateError, match="at least one required family"):
        S.Substrate(substrate_id="x", label="X", required_families=())


def test_a_substrate_with_no_discriminator_is_refused():
    with pytest.raises(S.SubstrateError, match="discriminates it from its class"):
        S.Substrate(substrate_id="x", label="X", required_families=("GH1",))


def test_declaring_sharing_without_a_caveat_is_refused():
    with pytest.raises(S.SubstrateError, match="what that sharing does not prove"):
        S.Substrate(
            substrate_id="x", label="X", required_families=("GH1",),
            discriminating_families=("GH2",), shares_with={"carb.inulin": "because"},
        )


# --------------------------------------------------------------------------- #
# the four named overclaims
# --------------------------------------------------------------------------- #


def test_gh13_alone_is_general_amylolysis_not_resistant_starch():
    """The specification's first named overclaim."""
    hits = _hits(("amyA", "GH13"))
    resistant = S.assess(S.BY_ID["carb.resistant_starch"], hits)
    assert resistant.specificity == "class_level_only"
    assert "nothing that distinguishes" in resistant.claim
    # ...while the general-starch card is entitled to say what it says.
    general = S.assess(S.BY_ID["carb.starch_general"], hits)
    assert general.specificity == "substrate_specific"


def test_resistant_starch_becomes_specific_only_with_its_own_subfamily():
    specific = S.assess(
        S.BY_ID["carb.resistant_starch"], _hits(("amyA", "GH13"), ("sas6", "GH13_36")),
    )
    assert specific.specificity == "substrate_specific"
    assert "GH13_36" in specific.discriminators_present


def test_gh32_alone_does_not_separate_fos_from_inulin():
    """The second named overclaim."""
    hits = _hits(("sacA", "GH32"))
    for substrate_id in ("carb.inulin", "carb.fos"):
        finding = S.assess(S.BY_ID[substrate_id], hits)
        assert finding.specificity == "class_level_only", substrate_id
    assert "does not separate inulin" in S.BY_ID["carb.inulin"].must_not_conclude


def test_gh2_and_gh42_alone_do_not_establish_gos_utilisation():
    """The third named overclaim."""
    finding = S.assess(S.BY_ID["carb.gos"], _hits(("lacZ", "GH2"), ("lacA", "GH42")))
    assert finding.specificity == "class_level_only"
    assert "Every lactose-using organism carries a beta-galactosidase" in \
        S.BY_ID["carb.gos"].must_not_conclude


def test_gh18_alone_does_not_demonstrate_chitin_utilisation():
    """The fourth named overclaim."""
    partial = S.assess(S.BY_ID["carb.chitin"], _hits(("chiA", "GH18")))
    assert partial.specificity == "class_level_only"
    assert "peptidoglycan" in S.BY_ID["carb.chitin"].must_not_conclude
    whole = S.assess(S.BY_ID["carb.chitin"], _hits(("chiA", "GH18"), ("nagZ", "GH20")))
    assert whole.specificity == "substrate_specific"


# --------------------------------------------------------------------------- #
# grading
# --------------------------------------------------------------------------- #


def test_nothing_found_is_insufficient_not_a_negative_result():
    finding = S.assess(S.BY_ID["carb.pectin"], [])
    assert finding.specificity == "insufficient"
    assert finding.required_missing == S.BY_ID["carb.pectin"].required_families
    assert finding.fragments == 0


def test_not_assayed_is_distinct_from_found_nothing():
    finding = S.assess(S.BY_ID["carb.pectin"], [], assayed=False)
    assert finding.specificity == "not_assayed"
    assert "was not searched for" in finding.claim


def test_a_subfamily_satisfies_its_parent_but_not_the_reverse():
    assert "GH13" in S.families_found([S.GeneHit("x", "GH13_36", 1)])
    assert "GH13_36" not in S.families_found([S.GeneHit("x", "GH13", 1)])


def test_every_finding_repeats_the_caveat_that_belongs_to_it():
    for substrate in S.SUBSTRATES:
        payload = S.assess(substrate, _hits(("g", substrate.required_families[0]))).to_json()
        limits = " ".join(payload["limitations"])
        assert "not a measured amount" in limits
        if substrate.must_not_conclude:
            assert substrate.must_not_conclude in limits


# --------------------------------------------------------------------------- #
# the arithmetic
# --------------------------------------------------------------------------- #


def test_a_shared_gene_is_counted_once_in_a_total():
    """The rule that makes shared evidence safe to show twice."""
    hits = [
        S.GeneHit("susB", "GH13", 100),   # counted for starch and for IMO
        S.GeneHit("susB", "GH31", 100),
        S.GeneHit("other", "GH13", 50),
    ]
    assert S.aggregate_fragments(hits) == 150
    assert sum(h.fragments for h in hits) == 250, "the naive sum is the error being avoided"


def test_the_shared_gene_report_names_the_arithmetic_risk():
    findings = S.assess_all({
        "GH13": [S.GeneHit("amyA", "GH13", 100)],
        "GH31": [S.GeneHit("amyA", "GH31", 100)],
    })
    report = S.shared_gene_report(findings)
    assert report["n_shared_genes"] >= 1
    assert "amyA" in report["shared"]
    assert set(report["shared"]["amyA"]) >= {"carb.starch_general", "carb.imo"}
    assert "counts each gene once" in report["note"]


def test_assess_all_produces_one_card_per_substrate():
    findings = S.assess_all({"GH32": [S.GeneHit("sacA", "GH32", 40)]})
    assert len(findings) == len(S.SUBSTRATES)
    by_id = {f.substrate.substrate_id: f for f in findings}
    assert by_id["carb.inulin"].specificity == "class_level_only"
    assert by_id["carb.pectin"].specificity == "insufficient"


def test_assess_all_routes_a_subfamily_hit_to_the_parent_card():
    findings = {f.substrate.substrate_id: f for f in S.assess_all({
        "GH13_36": [S.GeneHit("sas6", "GH13_36", 70)],
    })}
    assert findings["carb.resistant_starch"].specificity == "substrate_specific"
    assert findings["carb.starch_general"].specificity == "substrate_specific"


# --------------------------------------------------------------------------- #
# reconciliation with dbCAN
# --------------------------------------------------------------------------- #


def test_the_dbcan_mapping_is_read_by_column_name(tmp_path):
    path = tmp_path / "fam-substrate-mapping.tsv"
    path.write_text(
        "Substrate_high_level\tSubstrate_Simple\tFamily\tName\n"
        "cellulose\tcellulose\tGH5\tendoglucanase\n"
        "starch\tstarch\tGH13\talpha-amylase\n",
        encoding="utf-8",
    )
    mapping = S.load_dbcan_mapping(path)
    assert "cellulose" in mapping["GH5"]
    assert "starch" in mapping["GH13"]


def test_a_mapping_without_the_expected_columns_is_refused(tmp_path):
    path = tmp_path / "bad.tsv"
    path.write_text("a\tb\n1\t2\n", encoding="utf-8")
    with pytest.raises(S.SubstrateError, match="expected a family column"):
        S.load_dbcan_mapping(path)


def test_disagreement_with_dbcan_is_reported_not_resolved():
    mapping = {"GH5": frozenset({"cellulose"}), "GH13": frozenset({"starch"})}
    out = S.validate_against_dbcan(mapping)
    assert out["agrees"] is False
    assert out["families_not_in_mapping"], "families absent from this cut-down mapping"
    assert "question for a person" in out["note"]


def test_agreement_is_reported_when_the_mapping_covers_everything():
    mapping: dict[str, frozenset[str]] = {}
    for substrate in S.SUBSTRATES:
        words = S.DBCAN_SUBSTRATE_WORDS[substrate.substrate_id]
        for family in (*substrate.required_families, *substrate.discriminating_families):
            mapping.setdefault(family, frozenset())
            mapping[family] = mapping[family] | {words[0]}
    out = S.validate_against_dbcan(mapping)
    assert out["agrees"] is True
    assert out["families_not_in_mapping"] == {}


# --------------------------------------------------------------------------- #
# A02 — fermentation routes and the cross-feeding network
# --------------------------------------------------------------------------- #

from openbiota.extension import fermentation as F  # noqa: E402


def test_every_route_declares_what_it_produces_or_consumes():
    for route in F.ROUTES:
        assert route.produces or route.consumes, route.route_id
        assert route.required_genes, route.route_id


def test_a_route_without_required_genes_is_refused():
    with pytest.raises(F.FermentationError, match="needs required genes"):
        F.Route(route_id="x", label="X", required_genes=(), produces=("acetate",))


def test_a_route_that_neither_produces_nor_consumes_is_refused():
    with pytest.raises(F.FermentationError, match="what it produces or consumes"):
        F.Route(route_id="x", label="X", required_genes=("a",))


def test_a_gene_cannot_be_both_required_and_insufficient_alone():
    with pytest.raises(F.FermentationError, match="cannot be both"):
        F.Route(route_id="x", label="X", required_genes=("acs",),
                produces=("acetate",), insufficient_alone={"acs": "because"})


def test_acs_supports_assimilation_not_formation():
    """The specification names this one: an AMP-forming acs is not acetate formation."""
    finding = F.assess_route(F.BY_ID["acetate.pta_ack"], ["acs"])
    assert finding.state == "insufficient"
    assert "acs" in finding.insufficient_present
    payload = finding.to_json()
    assert "assimilation" in payload["insufficient_alone_present"]["acs"]
    assert F.assess_route(F.BY_ID["acetate.pta_ack"], ["pta"]).state == "supported"


def test_a_glycyl_radical_enzyme_is_not_isla():
    """The specification names this one too."""
    finding = F.assess_route(F.BY_ID["sulfur.taurine"], ["grdB"])
    assert finding.state == "insufficient"
    assert "activase" in finding.to_json()["insufficient_alone_present"]["grdB"]


def test_a_partly_present_route_is_partial_not_absent():
    # The succinate route needs its dehydrogenase and its reductase.
    finding = F.assess_route(F.BY_ID["succinate.to_propionate"], ["mmdA"])
    assert finding.state == "partial"
    assert finding.present == ("mmdA",) and finding.missing == ("scpA",)


def test_acetate_rests_on_the_committed_step_not_on_a_decoy_elsewhere():
    """ackA is the butyrate panel's decoy; duplicating it would take its reads."""
    route = F.BY_ID["acetate.pta_ack"]
    assert route.required_genes == ("pta",)
    assert "ackA" in route.corroborating_genes
    assert F.assess_route(route, ["pta"]).state == "supported"
    assert "would take reads from that panel" in route.note


def test_nothing_searched_is_distinct_from_nothing_found():
    assert F.assess_route(F.BY_ID["acetate.pta_ack"], [], assayed=False).state == "not_assayed"
    assert F.assess_route(F.BY_ID["acetate.pta_ack"], []).state == "absent"


def test_propionate_routes_stay_three_distinct_chemistries():
    ids = {r.route_id for r in F.ROUTES if "propionate" in r.produces}
    assert ids == {
        "lactate.utilisation_propionate", "succinate.to_propionate",
        "propionate.propanediol",
    }


def test_sulfur_routes_keep_their_separate_endpoints():
    sulfide = [r for r in F.ROUTES if "hydrogen_sulfide" in r.produces]
    assert len(sulfide) >= 2
    assert {r.route_id for r in sulfide} >= {"hydrogen.consumption_sulfate", "sulfur.taurine"}
    assert "not every sulfur pathway" in F.BY_ID["sulfur.taurine"].note


def test_production_and_consumption_are_never_netted():
    view = F.summarise(["hydA", "mcrA"])
    hydrogen = view["metabolites"]["hydrogen"]
    assert hydrogen["produced_by"] and hydrogen["consumed_by"]
    assert "not netted" in hydrogen["note"]
    assert "no universal adverse direction" in " ".join(view["limitations"])


def test_only_supported_routes_draw_a_solid_edge():
    view = F.summarise(["pta", "ackA", "hydA"])
    styles = {e["style"] for e in view["network"]["edges"]}
    assert "solid" in styles
    partial = F.summarise(["mmdA"])  # succinate-to-propionate incomplete
    propionate = [
        e for e in partial["network"]["edges"]
        if e["to"] == "propionate" and e["route_id"] == "succinate.to_propionate"
    ]
    assert propionate and all(e["style"] == "dashed" for e in propionate)


def test_the_network_disclaims_being_a_flow():
    limits = " ".join(F.summarise(["pta", "ackA"])["network"]["limitations"])
    assert "not flows" in limits
    assert "does not assume its substrate is present" in limits
    assert "would establish there" in limits


def test_the_existing_scores_are_named_as_preserved():
    view = F.summarise([])
    assert set(view["preserved_unchanged"]) >= {"butyrate", "propionate", "methane"}
    assert "does not rescore them" in view["preserved_note"]


# --------------------------------------------------------------------------- #
# A03 — the nine-vitamin dashboard
# --------------------------------------------------------------------------- #

from openbiota.extension import vitamins as V  # noqa: E402


def test_nine_vitamins_four_new_and_five_reused():
    assert len(V.NEW_VITAMINS) == 4
    assert len(V.REUSED) == 5
    assert {v.vitamin_id for v in V.NEW_VITAMINS} == {"b1", "b3", "b5", "b6"}
    assert {r["vitamin_id"] for r in V.REUSED} == {"b2", "b7", "b9", "b12", "k2"}


def test_a_transport_gene_is_never_counted_as_synthesis():
    """The distinction the three columns exist to keep."""
    b1 = V.NEW_VITAMINS[0]
    finding = V.assess_vitamin(b1, ["thiB", "thiP", "thiQ"])
    assert finding.by_column["uptake"] == "complete"
    assert finding.by_column["synthesis"] == "absent"
    assert finding.uptake_only and not finding.can_synthesise
    assert "can take it in, but none carries a route to make it" in finding.headline
    assert "never added together" in finding.to_json()["columns_note"]


def test_one_branch_is_not_a_pathway():
    b1 = V.NEW_VITAMINS[0]
    thiazole_only = V.assess_vitamin(b1, ["thiG", "thiH", "thiS"])
    assert thiazole_only.by_column["synthesis"] == "partial"
    both_branches = V.assess_vitamin(b1, ["thiG", "thiH", "thiS", "thiC", "thiD"])
    assert both_branches.by_column["synthesis"] == "partial", "uncoupled branches"
    coupled = V.assess_vitamin(b1, ["thiG", "thiH", "thiS", "thiC", "thiD", "thiE"])
    assert coupled.by_column["synthesis"] == "complete"


def test_a_multi_branch_synthesis_route_must_declare_its_coupling():
    with pytest.raises(V.VitaminError, match="must name its coupling"):
        V.VitaminRoute(
            route_id="x", label="X", column="synthesis",
            branches=(V.Branch("a", "A", ("g1",)), V.Branch("b", "B", ("g2",))),
        )


def test_b6_alternatives_are_alternatives_not_requirements():
    b6 = next(v for v in V.NEW_VITAMINS if v.vitamin_id == "b6")
    assert V.assess_vitamin(b6, ["pdxS", "pdxT"]).by_column["synthesis"] == "complete"
    assert V.assess_vitamin(b6, ["pdxA", "pdxJ"]).by_column["synthesis"] == "complete"
    # ...but one isolated homolog is not a route.
    assert V.assess_vitamin(b6, ["pdxS"]).by_column["synthesis"] == "partial"
    assert "alternatives" in b6.must_not_conclude


def test_downstream_use_is_not_a_second_production_route():
    b5 = next(v for v in V.NEW_VITAMINS if v.vitamin_id == "b5")
    genes = {g for r in b5.routes for g in r.all_genes}
    assert "coaA" not in genes, "CoA synthesis consumes pantothenate; it does not make it"
    assert "coenzyme A" in b5.must_not_conclude


def test_b3_keeps_de_novo_and_salvage_apart():
    b3 = next(v for v in V.NEW_VITAMINS if v.vitamin_id == "b3")
    salvage_only = V.assess_vitamin(b3, ["pncA", "pncB"])
    assert salvage_only.by_column["salvage"] == "complete"
    assert salvage_only.by_column["synthesis"] == "absent"
    assert not salvage_only.can_synthesise
    assert "rebuild it from a partial form" in salvage_only.headline


def test_every_vitamin_has_a_synthesis_route_to_compare_uptake_against():
    with pytest.raises(V.VitaminError, match="no synthesis route"):
        V.Vitamin(
            vitamin_id="x", label="X", common_name="x",
            routes=(V.VitaminRoute("r", "R", "uptake", (V.Branch("b", "B", ("g",)),)),),
        )


def test_nothing_searched_is_distinct_from_no_route_found():
    b1 = V.NEW_VITAMINS[0]
    assert V.assess_vitamin(b1, [], assayed=False).by_column["synthesis"] == "not_assayed"
    assert V.assess_vitamin(b1, []).by_column["synthesis"] == "absent"


def test_the_dashboard_reuses_the_five_without_rescoring_them():
    view = V.dashboard([], existing_panels={"riboflavin": {"copies_per_100_genomes": 12.3}})
    assert view["n_vitamins"] == 9
    riboflavin = next(r for r in view["reused"] if r["vitamin_id"] == "b2")
    assert riboflavin["reused_unchanged"] is True
    assert riboflavin["existing_value"] == 12.3
    assert "without rescoring" in riboflavin["note"]


def test_the_reused_route_details_name_their_own_traps():
    details = {r["vitamin_id"]: r["route_detail"] for r in V.REUSED}
    assert "more than one corrinoid" in details["b12"]
    assert "do not identify a chain length" in details["k2"]


def test_the_dashboard_disclaims_being_a_blood_test():
    limits = " ".join(V.dashboard([])["limitations"])
    assert "not serum vitamin status" in limits
    assert "not host absorption" in limits
    assert "stop a prescribed supplement" in limits


def test_the_gene_inventory_records_every_imported_gene_with_its_route():
    inventory = V.imported_gene_inventory()
    assert len(inventory) > 20
    genes = {row["gene"] for row in inventory}
    assert {"thiE", "panC", "nadA", "pdxS"} <= genes
    coupling = [r for r in inventory if r["role"] == "coupling"]
    assert {r["gene"] for r in coupling} >= {"thiE", "panC", "nadC"}
    for row in inventory:
        assert row["column"] in V.COLUMN_LABELS


def test_a_gene_no_panel_searches_for_is_not_assayed_not_absent():
    """A gap in the assay is not a property of the person."""
    b1 = V.NEW_VITAMINS[0]
    # The run searched for butyrate genes and found them; it never looked
    # for thiamine genes at all.
    finding = V.assess_vitamin(b1, ["but", "buk"], genes_searched=["but", "buk", "ackA"])
    assert finding.by_column["synthesis"] == "not_assayed"
    assert not finding.searched
    assert "were not searched for" in finding.headline
    # Searched and genuinely not found is a different statement.
    looked = V.assess_vitamin(b1, [], genes_searched=["thiG", "thiH", "thiS", "thiC"])
    assert looked.by_column["synthesis"] == "absent"
    assert looked.searched


def test_a_column_the_vitamin_has_no_route_for_is_not_applicable():
    b5 = next(v for v in V.NEW_VITAMINS if v.vitamin_id == "b5")
    assert not b5.routes_in("salvage")
    finding = V.assess_vitamin(b5, [], genes_searched=["panB"])
    assert finding.by_column["salvage"] == "not_applicable", (
        "no modelled salvage route is not the same as an absent one"
    )
    assert "not_applicable" in V.ROUTE_STATES


# --------------------------------------------------------------------------- #
# A04 — protein, nitrogen and aromatic metabolism
# --------------------------------------------------------------------------- #

from openbiota.extension import nitrogen as N  # noqa: E402


def test_all_six_nitrogen_modules_are_present():
    assert len(N.MODULE_DEFINITIONS) == 6
    assert {m.module_id for m in N.MODULE_DEFINITIONS} == set(N.MODULE_LABELS)


def test_protein_breakdown_bcaa_synthesis_and_bcaa_fermentation_stay_three_things():
    """The specification names this confusion explicitly."""
    view = N.summarise([])
    concepts = view["three_different_concepts"]
    assert {"protein_breakdown", "bcaa_synthesis", "bcaa_fermentation"} <= set(concepts)
    assert "three different things" in concepts["note"]
    bcfa = N.BY_ID["bcfa"]
    assert "opposite directions" in bcfa.distinct_from["BCAA biosynthesis"]
    assert "bcaa" in bcfa.links_to, "it links to the existing panel"
    assert bcfa.module_id != "bcaa", "and does not replace it"


def test_generic_peptidase_abundance_is_not_protein_digestion():
    proteolysis = N.BY_ID["proteolysis"]
    assert not proteolysis.aggregatable
    assert "not how much of your" in proteolysis.distinct_from["dietary protein digestion"]
    payload = N.assess_module(proteolysis, [], searched=list(proteolysis.genes)).to_json()
    assert "would dominate the total" in payload["no_aggregate_reason"]


def test_housekeeping_genes_do_not_carry_the_finding():
    proteolysis = N.BY_ID["proteolysis"]
    housekeeping_only = N.assess_module(
        proteolysis, ["clpP", "lon", "pepA", "pepN"], searched=list(proteolysis.genes),
    )
    assert housekeeping_only.informative_present == ()
    assert housekeeping_only.state == "partial", "not 'present' on housekeeping alone"


def test_a_housekeeping_gene_must_belong_to_its_module():
    with pytest.raises(N.NitrogenError, match="marked housekeeping but not in"):
        N.NitrogenModule(
            module_id="ammonia", label="X", genes=("gdhA",), housekeeping=("nope",),
        )


def test_a_microbial_precursor_is_never_the_host_metabolite():
    steps = {h.microbial_product: h for h in N.HOST_STEPS}
    pagln = steps["phenylacetate"]
    assert pagln.host_product == "phenylacetylglutamine (PAGln)"
    assert pagln.to_json()["not_measured_here"] == "phenylacetylglutamine (PAGln)"
    assert "your liver conjugates it" in pagln.note
    tmao = steps["trimethylamine (TMA)"]
    assert "FMO3" in tmao.host_enzyme
    assert "varies widely between people" in tmao.note


def test_urease_stays_its_own_reading():
    ammonia = N.BY_ID["ammonia"]
    assert "urease" in ammonia.links_to
    assert "urease" not in ammonia.genes
    assert "unchanged" in ammonia.distinct_from["urease"]
    assert "urease" in N.summarise([])["preserved_unchanged"]


def test_the_aromatic_module_links_without_merging_the_existing_readings():
    aromatic = N.BY_ID["aromatic"]
    assert set(aromatic.links_to) == {"indole", "ipa", "pcresol"}
    assert "stay separate" in aromatic.note
    assert aromatic.host_steps


def test_a_module_whose_genes_were_never_searched_is_not_assayed():
    view = N.summarise(["ilvE", "kivD"], searched=["ilvE", "kivD", "bkdA"])
    states = {m["module_id"]: m["state"] for m in view["modules"]}
    assert states["bcfa"] == "partial"
    assert states["ammonia"] == "not_assayed", "no ammonia gene was searched for"


# --------------------------------------------------------------------------- #
# A05, A06, A08 — neuroactive, plant-compound and surface biotransformations
# --------------------------------------------------------------------------- #

from openbiota.extension import biotransform as B  # noqa: E402


def test_every_step_says_what_it_does_not_establish():
    for step in B.ALL_STEPS:
        assert step.does_not_establish, step.step_id
        assert step.feature_id in B.FEATURES
    with pytest.raises(B.BiotransformError, match="does not establish"):
        B.Step(step_id="x", feature_id="A05", label="X",
               requirement=B.Requirement(all_of=("g",)),
               establishes="something", does_not_establish="")


def test_the_glucosinolate_core_is_the_boolean_the_spec_states():
    """BT2158 AND (BT2156 OR BT2157), not five genes flattened together."""
    step = B.BY_ID["diet.glucosinolate_isothiocyanate_conversion"]
    assert step.requirement.describe() == "BT2158 and (BT2156 or BT2157)"
    satisfied = step.requirement.satisfied_by
    assert satisfied(frozenset({"BT2158", "BT2156"}))
    assert satisfied(frozenset({"BT2158", "BT2157"}))
    assert not satisfied(frozenset({"BT2158"})), "the OR arm is required"
    assert not satisfied(frozenset({"BT2156", "BT2157"})), "BT2158 is required"


def test_a_requirement_must_require_something():
    with pytest.raises(B.BiotransformError, match="must require something"):
        B.Requirement()


def test_a_matched_negative_control_withholds_the_step():
    """UrdA is close enough to be mistaken for a urolithin dehydroxylase."""
    step = B.BY_ID["polyphenol.urolithin_9_dehydroxylation"]
    assert "urdA" in step.negative_controls
    clean = B.assess_step(step, ["ucdC", "ucdF", "ucdO"])
    assert clean.state == "supported"
    contaminated = B.assess_step(step, ["ucdC", "ucdF", "ucdO", "urdA"])
    assert contaminated.state == "negative_control_matched"
    payload = contaminated.to_json()
    assert "urocanate" in payload["negative_controls_matched"]["urdA"]
    assert "signal about the panel, not a finding about you" in " ".join(payload["limitations"])


def test_a_gene_cannot_be_both_required_and_a_negative_control():
    with pytest.raises(B.BiotransformError, match="both required and a negative control"):
        B.Step(step_id="x", feature_id="A06", label="X",
               requirement=B.Requirement(all_of=("ucdC",)),
               establishes="a", does_not_establish="b",
               negative_controls={"ucdC": "why"})


def test_all_three_urolithin_components_are_required():
    step = B.BY_ID["polyphenol.urolithin_9_dehydroxylation"]
    assert set(step.requirement.all_of) == {"ucdC", "ucdF", "ucdO"}
    assert B.assess_step(step, ["ucdC", "ucdF"]).state == "partial"
    assert "whole ellagitannin-to-urolithin pathway" in step.does_not_establish


def test_dna_evidence_is_never_a_metabotype():
    for step_id in ("polyphenol.urolithin_9_dehydroxylation",
                    "polyphenol.equol_daidzein_conversion"):
        assert "metabotype" in B.BY_ID[step_id].does_not_establish
    limits = " ".join(B.summarise([])["limitations"])
    assert "A DNA finding is not a metabotype" in limits
    assert "substrate challenge" in limits


def test_a_general_reductase_domain_is_not_equol_chemistry():
    step = B.BY_ID["polyphenol.equol_daidzein_conversion"]
    assert set(step.requirement.all_of) == {"dznr", "ddr", "tdr"}
    assert "generic_reductase" in step.negative_controls
    assert "eqlA" in step.negative_controls["generic_reductase"]


def test_gaba_production_and_degradation_are_never_subtracted():
    balance = B.gaba_balance(10.0, 4.0)
    assert balance["synthesis"] == 10.0 and balance["degradation"] == 4.0
    assert balance["ratio"] == pytest.approx(0.4)
    assert balance["formula"].startswith("degradation / synthesis")
    assert "is not a net flux" in balance["note"]
    assert "difference" in balance["note"]


def test_a_ratio_is_withheld_rather_than_dividing_by_zero():
    balance = B.gaba_balance(0.0, 5.0)
    assert balance["ratio"] is None
    assert "divide by zero" in balance["ratio_withheld"]


def test_mucin_capacity_is_not_barrier_damage():
    step = B.BY_ID["mucin.glycan_foraging"]
    assert "not an image of" in step.does_not_establish
    assert "physiological" in step.does_not_establish
    assert "dietary-fibre CAZymes" in step.note


def test_lipid_a_genes_do_not_resolve_a_tlr4_stimulus():
    step = B.BY_ID["lps.lipid_a_modification"]
    assert "TLR4" in step.does_not_establish
    assert "total Gram-negative abundance resolves nothing" in step.does_not_establish
    assert "fabricated endotoxin" in step.note


def test_bile_chemistry_is_a_graph_not_a_chain():
    hsdh = B.BY_ID["bile.hsdh_transformation"]
    assert "graph rather than a chain" in hsdh.does_not_establish
    bsh = B.BY_ID["bile.deconjugation"]
    assert "acyltransferase" in bsh.does_not_establish, "BSH can also reconjugate"
    assert "rumen validation population supplies no" in hsdh.note


def test_a_step_whose_genes_were_never_searched_is_not_assayed():
    step = B.BY_ID["bile.dehydroxylation"]
    assert B.assess_step(step, [], searched=["gadA"]).state == "not_assayed"
    assert B.assess_step(step, [], searched=["baiB"]).state == "absent"


# --------------------------------------------------------------------------- #
# A07 — GLP-1-related microbial mechanisms
# --------------------------------------------------------------------------- #


def test_the_four_glp1_channels_stay_four_separate_mechanisms():
    view = B.glp1_panel([])
    assert view["n_channels"] == 4
    assert len({c["step_id"] for c in view["channels"]}) == 4
    assert "not summed into a single" in " ".join(view["limitations"])


def test_host_receptors_are_named_as_not_microbial():
    """A stool metagenome cannot contain FFAR2."""
    targets = B.glp1_panel([])["not_microbial_targets"]
    assert {"FFAR2", "FFAR3", "GCG", "TGR5", "DPP4"} <= set(targets)
    assert "your microbes do not" in targets["FFAR2"]
    assert "human proglucagon" in targets["GCG"]
    assert "different proteins with different substrates" in targets["DPP4"]
    # ...and none of them is a gene any channel searches for.
    searched = {g for step in B.GLP1_CHANNELS for g in step.requirement.genes}
    assert not (searched & set(targets))


def test_the_host_boundary_is_stated_not_implied():
    boundary = B.glp1_panel([])["host_boundary"]
    assert "human hormone" in boundary
    assert "nothing in this panel is a GLP-1 measurement" in boundary
    assert "A stool metagenome contains none of them" in boundary


def test_indole_routes_carry_no_universal_direction():
    """Section 9.1 says so explicitly."""
    indole = next(s for s in B.GLP1_CHANNELS if s.step_id == "glp1.indole_signalling")
    assert "no universal" in indole.does_not_establish
    assert "both stimulation and inhibition" in indole.does_not_establish
    assert "no universal polarity" in indole.note


def test_the_glp1_panel_reads_the_existing_scfa_readings_without_changing_them():
    scfa = next(s for s in B.GLP1_CHANNELS if s.step_id == "glp1.scfa_signalling")
    assert set(scfa.links_to) == {"butyrate", "propionate"}
    assert "unchanged" in scfa.note


# --------------------------------------------------------------------------- #
# the shared database: a duplicate target is a subtraction
# --------------------------------------------------------------------------- #


def test_no_gene_symbol_is_claimed_by_two_panels():
    """A duplicate target does not add a reading; it takes the references.

    This cost the IPA panel all but 44 of its 1,297 references, dropped
    butyrate by 5%, and took the hydrogen-sulfide panel's dsrA to zero
    fragments, which flipped a concordance call. Alignment is competitive
    across one database, so two panels cannot both own a gene.
    """
    import collections
    from pathlib import Path

    import yaml

    panels_dir = Path(__file__).resolve().parents[1] / "panels"
    owner: dict[str, list[str]] = collections.defaultdict(list)
    for path in sorted(panels_dir.glob("*.yaml")):
        if path.stem == "_normalizer":
            continue
        spec = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        for kind in ("targets", "decoys"):
            for entry in spec.get(kind) or []:
                gene = str(entry.get("gene") or "").strip().lower()
                if gene:
                    owner[gene].append(f"{path.stem}:{kind}:{entry['id']}")
    clashes = {gene: rows for gene, rows in owner.items() if len(rows) > 1}
    assert not clashes, (
        "these gene symbols are claimed by more than one panel, so they compete "
        f"for the same reads: {clashes}"
    )


def test_every_panel_target_declares_a_gene_so_the_clash_check_can_see_it():
    from pathlib import Path

    import yaml

    panels_dir = Path(__file__).resolve().parents[1] / "panels"
    missing: list[str] = []
    for path in sorted(panels_dir.glob("*.yaml")):
        if path.stem == "_normalizer":
            continue
        spec = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
        for entry in spec.get("targets") or []:
            if not str(entry.get("gene") or "").strip():
                missing.append(f"{path.stem}:{entry['id']}")
    assert not missing, f"targets without a gene symbol escape the clash check: {missing}"
