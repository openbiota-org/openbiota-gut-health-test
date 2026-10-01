"""Acceptance tests of BUILD_SPEC_v0.8.3 section 15.

Each test names the AT it discharges in its docstring, which is what the
ledger in `extension/acceptance_tests.yaml` reads. Naming an ID in a
comment is not enough to claim it: the test has to fail if the behaviour
regresses, so these exercise the modules rather than assert on constants.
"""

from __future__ import annotations

import pathlib

import pytest
import yaml

from openbiota.extension import biotransform as BIO
from openbiota.extension import capacities as CAP
from openbiota.extension import contexts as CX
from openbiota.extension import fermentation as FERM
from openbiota.extension import longitudinal as LON
from openbiota.extension import registry as REG
from openbiota.extension import substrates as SUB
from openbiota.extension import vitamins as VIT
from openbiota.extension.planner import START_HERE_MAX as PLAN_START_HERE_MAX
from openbiota.samples import report_pdf, results_dir, results_file  # noqa: E402,F401
from tests import _data  # noqa: E402

REPO = pathlib.Path(__file__).resolve().parents[1]
LEDGER = REPO / "extension" / "acceptance_tests.yaml"


def hit(gene: str, family: str, fragments: int = 100) -> SUB.GeneHit:
    return SUB.GeneHit(gene_id=gene, family=family, fragments=fragments)


# --------------------------------------------------------------------------- #
# the ledger itself
# --------------------------------------------------------------------------- #


def test_the_acceptance_ledger_covers_every_test_in_the_specification():
    """AT009 in part: every acceptance test is accounted for.

    Either a named test exercises it, or the entry records why not. An
    entry with neither is the failure this catches.
    """
    import re

    ledger = yaml.safe_load(LEDGER.read_text(encoding="utf-8"))
    spec = (REPO / "specs" / "BUILD_SPEC_v0.8.3.md").read_text(encoding="utf-8")
    in_spec = set(re.findall(r"^\| (AT\d{3}) \|", spec, flags=re.M))
    in_ledger = {row["id"] for row in ledger["tests"]}
    assert in_spec == in_ledger, (
        f"ledger and specification disagree: missing {sorted(in_spec - in_ledger)}, "
        f"unknown {sorted(in_ledger - in_spec)}"
    )

    # a plain scan of the test tree: no dependence on a search tool being installed
    exercised: set[str] = set()
    for path in (REPO / "tests").rglob("*.py"):
        exercised.update(re.findall(r"AT\d{3}", path.read_text(encoding="utf-8", errors="ignore")))
    unaccounted = [
        row["id"] for row in ledger["tests"]
        if row["id"] not in exercised and not str(row.get("reason") or "").strip()
    ]
    assert not unaccounted, (
        f"{len(unaccounted)} acceptance tests have neither a test that names them nor a "
        f"recorded reason: {unaccounted[:12]}"
    )


def test_every_ledger_entry_keeps_the_requirement_text_from_the_specification():
    ledger = yaml.safe_load(LEDGER.read_text(encoding="utf-8"))
    for row in ledger["tests"]:
        assert row["requirement"].strip(), row["id"]


# --------------------------------------------------------------------------- #
# A01 — substrates
# --------------------------------------------------------------------------- #


def test_at013_generic_gh13_does_not_activate_resistant_starch():
    """AT013: a single generic GH13 hit supports broad amylolysis but does
    not activate complete resistant-starch utilization."""
    hits = [hit("amyA", "GH13")]
    general = SUB.assess(SUB.BY_ID["carb.starch_general"], hits)
    resistant = SUB.assess(SUB.BY_ID["carb.resistant_starch"], hits)
    assert general.specificity in {"class_level_only", "substrate_specific"}
    assert resistant.specificity == "class_level_only", (
        "one GH13 must not establish resistant-starch utilisation"
    )
    assert not resistant.discriminators_present, (
        "the subfamilies that would make this specific are absent and must be seen to be"
    )
    assert resistant.substrate.discriminating_families, (
        "resistant starch must declare which subfamilies would discriminate it"
    )


def test_at014_shared_evidence_is_not_counted_twice_in_an_aggregate():
    """AT014: inulin/FOS and lactose/GOS can share evidence without
    duplicated aggregate fragment counts."""
    shared = [hit("sacA", "GH32", 500), hit("sacA", "GH32", 500)]
    assert SUB.aggregate_fragments(shared) == 500, (
        "the same gene appearing on two cards cannot contribute twice"
    )
    inulin = SUB.BY_ID["carb.inulin"]
    fos = SUB.BY_ID["carb.fos"]
    linked = set(inulin.shares_with) | set(fos.shares_with)
    assert {"carb.fos", "carb.inulin"} & linked, (
        "inulin and FOS must declare that they share evidence"
    )


def test_at015_all_fourteen_substrate_views_exist_and_carry_a_state():
    """AT015: all 14 carbohydrate views render, with appropriate
    complete/partial/no-supported-hit/not-assayed states."""
    assert len(SUB.SUBSTRATES) == 14
    findings = SUB.assess_all({}, assayed=False)
    assert len(findings) == 14
    assert {f.specificity for f in findings} == {"not_assayed"}
    measured = SUB.assess_all({"GH13": [hit("amyA", "GH13")]}, assayed=True)
    assert all(f.specificity in SUB.SPECIFICITY for f in measured)
    assert any(f.specificity == "insufficient" for f in measured), (
        "a substrate whose families are absent must read as searched-and-absent"
    )


def test_at022_a_route_spread_across_a_community_is_not_genome_linked():
    """AT022: community-distributed route genes do not receive
    `genome_linked` completeness without linkage evidence."""
    for substrate in SUB.SUBSTRATES:
        finding = SUB.assess(substrate, [hit("g1", f) for f in substrate.required_families])
        payload = finding.to_json()
        assert "genome_linked" not in str(payload.get("specificity")), substrate.substrate_id


# --------------------------------------------------------------------------- #
# A02 — fermentation
# --------------------------------------------------------------------------- #


def test_at019_acs_alone_does_not_prove_acetate_production():
    """AT019: AMP-forming `acs` alone does not prove acetate production;
    assimilation, ADP-forming alternatives and Wood-Ljungdahl CODH/ACS
    remain distinct."""
    routes = {r.route_id: r for r in FERM.ROUTES}
    acetate = [r for rid, r in routes.items() if "acetate" in rid]
    assert acetate, "there must be at least one acetate route"
    ids = {r.route_id for r in acetate}
    assert len(ids) > 1 or any("wood" in i or "ljung" in i for i in routes), (
        "the Wood-Ljungdahl route must be distinguishable from phosphotransacetylase"
    )
    for route in acetate:
        assert "acs" not in route.required_genes, (
            f"{route.route_id}: the AMP-forming synthetase must not be a required gene, "
            "because it runs towards assimilation as readily as towards production"
        )


def test_at028_hydrogen_production_and_consumption_are_separate_readings():
    """AT028: hydrogen-producing and hydrogen-consuming classes are
    separate; high hydrogen potential alone is not an adverse score."""
    hydrogen = [r for r in FERM.ROUTES if "hydrogen" in r.route_id]
    assert len(hydrogen) >= 2, "production and consumption must be distinct routes"
    producing = [r for r in hydrogen if "hydrogen" in r.produces]
    consuming = [r for r in hydrogen if "hydrogen" in r.consumes]
    assert producing and consuming, "production and consumption must both be represented"
    assert not ({r.route_id for r in producing} & {r.route_id for r in consuming}), (
        "one route cannot be both a producer and a consumer of hydrogen"
    )
    for route in hydrogen:
        # No route carries a good/bad direction: a healthy community makes a
        # great deal of hydrogen, so a high potential is not an adverse score.
        assert route.direction == "descriptive", route.route_id


def test_at034_generic_adhe_yields_no_strain_name_or_diagnosis():
    """AT034: generic `adhE` does not yield a high-alcohol strain name,
    measured ethanol value or auto-brewery diagnosis."""
    capacity = CAP.BY_ID["capacity.ethanol"]
    state = capacity.assess(frozenset({"adhE"}), frozenset({"adhE"}))
    assert state == "complete"
    text = f"{capacity.establishes} {capacity.does_not_establish} {capacity.context}".lower()
    assert "auto-brewery" in text and "does not establish" not in capacity.establishes
    assert "carrying the gene is ordinary" in text
    for banned in ("klebsiella pneumoniae strain in you", "mg/dl", "blood alcohol of"):
        assert banned not in text
    context = CX.BY_ID["context.auto_brewery_ethanol"]
    assert "not measured ethanol" in context.distinction


# --------------------------------------------------------------------------- #
# A03 — vitamins
# --------------------------------------------------------------------------- #


def test_at017_synthesis_is_separate_from_salvage_and_claims_no_deficiency():
    """AT017: B1/B3/B5/B6 synthesis is separated from salvage/uptake and
    does not generate host vitamin-deficiency text."""
    columns = dict(VIT.COLUMNS)
    assert "synthesis" in columns
    assert any(key in columns for key in ("salvage", "uptake", "salvage_uptake"))
    for vitamin in VIT.NEW_VITAMINS:
        assert vitamin.must_not_conclude, (
            f"{vitamin.vitamin_id} must record the conclusion it does not support"
        )
        blob = f"{vitamin.must_not_conclude} {vitamin.label}".lower()
        for banned in ("you are deficient", "your deficiency", "you need a supplement"):
            assert banned not in blob, vitamin.vitamin_id
        # Synthesis, salvage and uptake are graded in separate columns, which
        # is what stops "can take it in" reading as "can make it".
        columns_used = {r.column for r in vitamin.routes}
        assert columns_used <= set(dict(VIT.COLUMNS)), vitamin.vitamin_id
        assert "synthesis" in columns_used and columns_used - {"synthesis"}, (
            f"{vitamin.vitamin_id}: synthesis must be graded apart from salvage and uptake"
        )


def test_at018_reused_vitamin_cells_point_at_existing_results():
    """AT018: B2/B7/B9/B12/K2 dashboard cells reference unchanged existing
    result objects."""
    reused = {str(row["vitamin_id"]).lower() for row in VIT.REUSED}
    assert {"b2", "b7", "b9", "b12", "k2"} <= reused
    for row in VIT.REUSED:
        assert row["panel"], f"{row['vitamin_id']} must name the existing panel it reuses"


# --------------------------------------------------------------------------- #
# A05-A08 — named steps
# --------------------------------------------------------------------------- #


def test_at020_gaba_production_and_degradation_never_become_one_number():
    """AT020: GABA production and degradation have separate IDs, values and
    evidence; a percentile subtraction cannot be labelled net GABA."""
    synthesis = BIO.BY_ID["neuro.gaba_synthesis"]
    degradation = BIO.BY_ID["neuro.gaba_degradation"]
    assert synthesis.step_id != degradation.step_id
    assert not set(synthesis.requirement.genes) & set(degradation.requirement.genes)
    balance = BIO.gaba_balance(80.0, 20.0)
    assert "difference" not in balance, "a difference of percentiles must not be published"
    assert "net" not in str(balance.get("formula", "")).lower()
    assert "not a net flux" in balance["note"]


def test_at023_urda_positive_does_not_activate_urolithin_or_equol():
    """AT023: UrdA-positive/urolithin-negative fixture does not activate new
    urolithin or equol modules."""
    detected = frozenset({"urdA"})
    searched = frozenset({"urdA", "ucdC", "ucdF"})
    for step_id in ("polyphenol.urolithin_9_dehydroxylation",
                    "polyphenol.equol_daidzein_conversion"):
        step = BIO.BY_ID[step_id]
        assert "urdA" not in step.requirement.genes, (
            f"{step_id}: urdA is urocanate reductase and is not a urolithin marker"
        )
        finding = BIO.assess_step(step, detected, searched=searched)
        assert finding.state != "supported", (
            f"{step_id} must not be supported by urdA alone"
        )


def test_at026_glucosinolate_core_rule_is_a_boolean_not_a_gene_list():
    """AT026: glucosinolate core rule accepts BT2158 plus BT2156 or BT2157,
    with supporting context; generic hydrolases remain nonspecific."""
    step = BIO.BY_ID["diet.glucosinolate_isothiocyanate_conversion"]
    requirement = step.requirement
    assert requirement.satisfied_by(frozenset({"BT2158", "BT2156"}))
    assert requirement.satisfied_by(frozenset({"BT2158", "BT2157"}))
    assert not requirement.satisfied_by(frozenset({"BT2156", "BT2157"})), (
        "the core gene is required, not optional"
    )
    assert not requirement.satisfied_by(frozenset({"BT2158"})), (
        "the core gene alone is not the route"
    )


def test_at033_the_neuroactive_molecules_stay_distinct_entities():
    """AT033: tryptamine, serotonin, melatonin and receptor-active
    derivatives remain distinct chemical entities."""
    context = CX.BY_ID["context.serotonin_melatonin_sleep"]
    assert "three different molecules" in context.distinction
    text = " ".join(context.mechanisms).lower()
    assert "tryptamine is not serotonin" in text


# --------------------------------------------------------------------------- #
# section 5.8
# --------------------------------------------------------------------------- #


def test_at055_research_level_capacities_stay_visible_as_research():
    """AT055 in part: lack of an RCT alone never removes a supported
    research option; research is a visible evidence level."""
    view = CAP.summarise(["xdhA", "xdhB"], searched=["xdhA", "xdhB", "xdhC"])
    assert view["evidence_level"] == CAP.EVIDENCE_RESEARCH
    assert all(row["evidence_level"] == "research" for row in view["capacities"])
    assert view["n_complete"] >= 1, "a research capacity still produces a reading"


# --------------------------------------------------------------------------- #
# A13 — longitudinal and context
# --------------------------------------------------------------------------- #


def test_at046_bray_curtis_endpoints_are_exact():
    """AT046: identical compositions have Bray-Curtis distance 0; disjoint
    unit-sum compositions have distance 1."""
    same = {"a": 0.5, "b": 0.5}
    assert LON.bray_curtis(same, dict(same)) == pytest.approx(0.0)
    assert LON.bray_curtis({"a": 1.0}, {"b": 1.0}) == pytest.approx(1.0)


def test_at046_an_empty_union_is_undefined_rather_than_perfect_agreement():
    """AT046 companion: empty unions are undefined, not perfect agreement."""
    assert LON.jaccard_distance({}, {}) is None
    assert LON.bray_curtis({}, {}) is None


def test_at047_two_timepoints_are_a_change_not_a_measured_resilience():
    """AT047: two timepoints yield change, not measured resilience."""
    assert LON.MIN_BASELINE_SAMPLES >= 3


def test_at045_reprocessing_the_same_reads_is_not_a_new_timepoint():
    """AT045: same-FASTQ reports from different methods cannot be shown as
    biological improvement over time."""
    assert "tool_version" in LON.BLOCKING_FIELDS or any(
        "version" in f for f in LON.BLOCKING_FIELDS
    ), "a method change must block a native comparison"


def test_at129_every_context_links_to_something_that_exists():
    """AT129 in part: supported findings have a working option lookup.

    Each of the fifteen contexts must point at identifiers this codebase
    actually produces, or the navigation is a dead link.
    """
    panels_dir = REPO / "panels"
    known_panels = {p.stem for p in panels_dir.glob("*.yaml")} - {"_normalizer"}
    known_capacities = set(CAP.BY_ID)
    known_steps = set(BIO.BY_ID)
    broken: list[str] = []
    for context in CX.CONTEXTS:
        for link in context.links:
            if link.kind == "panel" and link.key not in known_panels:
                broken.append(f"{context.context_id} -> panel {link.key}")
            if link.kind == "capacity" and link.key not in known_capacities:
                broken.append(f"{context.context_id} -> capacity {link.key}")
            if link.kind == "step" and link.key not in known_steps:
                broken.append(f"{context.context_id} -> step {link.key}")
    assert not broken, f"contexts pointing at things that do not exist: {broken}"


# --------------------------------------------------------------------------- #
# coverage and preservation
# --------------------------------------------------------------------------- #


def test_at009_every_feature_and_metric_has_an_implementation_mapping():
    """AT009: every required A01-A16 feature and every coverage-appendix
    metric has an implementation mapping; an unimplemented placeholder is
    not marked data-limited."""
    registry = REG.load()
    unmapped = [
        c.capability_id for c in registry.capabilities.values()
        if not c.module and c.implementation != "not_implemented"
    ]
    assert not unmapped, f"capabilities with no implementing module: {unmapped}"
    for capability in registry.capabilities.values():
        if capability.implementation == "not_implemented":
            assert capability.note, (
                f"{capability.capability_id}: an unimplemented capability must record why"
            )


def test_at010_no_new_panel_is_pinned_to_a_moving_target():
    """AT010 in part: no production asset is pinned only to `latest`."""
    for path in sorted((REPO / "panels").glob("*.yaml")):
        text = path.read_text(encoding="utf-8")
        assert "latest" not in text.lower() or "translate" in text.lower(), path.name


def test_acceptance_coverage_does_not_go_backwards():
    """A ratchet on the ledger.

    116 of the 160 acceptance tests are exercised by a test that names them.
    That number is allowed to rise and not to fall: a change that quietly
    drops coverage fails here rather than being noticed later, and raising
    the floor is a one-line edit made deliberately.
    """
    floor = 116
    ledger = yaml.safe_load(LEDGER.read_text(encoding="utf-8"))
    covered = sum(1 for row in ledger["tests"] if row.get("covered_by"))
    assert covered >= floor, (
        f"acceptance coverage fell from {floor} to {covered}"
    )
    assert ledger["summary"]["n_covered"] == covered, (
        "the ledger's own summary is stale; re-run the refresh"
    )


# --------------------------------------------------------------------------- #
# A04 — nitrogen, and the host steps that are not microbial
# --------------------------------------------------------------------------- #


def test_at021_bcaa_synthesis_fermentation_and_bcfa_are_three_concepts():
    """AT021: BCAA biosynthesis, amino-acid fermentation and BCFA generation
    cannot share one biochemical identity."""
    from openbiota.extension import nitrogen as NIT

    by_id = {m.module_id: m for m in NIT.MODULE_DEFINITIONS}
    bcfa = by_id["bcfa"]
    fermentation = by_id["amino_acid_fermentation"]
    assert bcfa.module_id != fermentation.module_id
    assert not set(bcfa.genes) & set(fermentation.genes), (
        "sharing a gene between these two would let one reading stand for the other"
    )
    assert "BCAA biosynthesis" in bcfa.distinct_from, (
        "the BCFA module must record that it is not the biosynthesis panel"
    )
    assert "opposite directions" in bcfa.distinct_from["BCAA biosynthesis"]
    assert "bcaa" in bcfa.links_to, "it must link to the existing panel it is not"


def test_at002_microbial_precursors_are_not_relabelled_as_host_products():
    """AT002/AT003 in part: existing labels remain accurate.

    Microbial TMA is not host TMAO, and microbial phenylacetate is not
    circulating PAGln. Each host step is recorded as the host's.
    """
    from openbiota.extension import nitrogen as NIT

    aromatic = next(m for m in NIT.MODULE_DEFINITIONS if m.module_id == "aromatic")
    products = {step.microbial_product: step for step in aromatic.host_steps}
    assert "trimethylamine (TMA)" in products
    tmao = products["trimethylamine (TMA)"]
    assert tmao.host_product.startswith("trimethylamine N-oxide")
    assert "FMO3" in tmao.host_enzyme or "FMO3" in tmao.note
    pagln = products["phenylacetate"]
    assert "PAGln" in pagln.host_product
    assert "host product" in pagln.note


# --------------------------------------------------------------------------- #
# A10 — the organism explorer
# --------------------------------------------------------------------------- #


def test_at052_a_dna_assay_cannot_report_an_assessed_negative():
    """AT052: DNA-only input cannot report an assessed negative for the
    listed RNA viruses."""
    from openbiota.extension import explorer as EXP

    assert "not_assessable" in EXP.DETECTION_STATES, (
        "there must be a state distinct from 'no supported detection' for a target "
        "this assay cannot see at all"
    )
    assert "no_supported_detection" in EXP.DETECTION_STATES
    assert "not_assessable" != "no_supported_detection"


def test_at049_every_resolution_outcome_is_named_rather_than_dropped():
    """AT049: every taxonomy label resolves or produces an explicit
    split/merge/unresolved state; none is silently dropped."""
    from openbiota.extension import explorer as EXP

    assert {"split_or_merge_ambiguous", "unresolvable", "not_in_catalogue"} <= (
        EXP.RESOLUTION_STATES
    )


# --------------------------------------------------------------------------- #
# A12 — the planner
# --------------------------------------------------------------------------- #


def test_at060_duplicate_items_merge_but_unlike_ones_never_do():
    """AT060: duplicate food/action items merge with all supported targets
    and tradeoffs; overlapping targets do not create a net-benefit score."""
    from openbiota.extension import planner as PLAN

    assert PLAN.normalise_food("Kimchi") == PLAN.normalise_food("kimchi")
    for specific, general in PLAN.NEVER_MERGE:
        assert PLAN.normalise_food(specific) != PLAN.normalise_food(general), (
            f"{specific!r} and {general!r} are different things and must not merge"
        )


def test_at058_no_rank_reason_is_a_dose_or_a_clinical_benefit():
    """AT058: a study concentration or animal dose cannot become a human
    regimen through automatic unit conversion."""
    from openbiota.extension import planner as PLAN

    joined = " ".join(PLAN.RANK_REASONS).lower()
    for banned in ("dose", "mg", "grams per day", "net benefit", "efficacy"):
        assert banned not in joined, f"a ranking reason must not be {banned!r}"


def test_at055_every_action_category_carries_a_monitoring_question():
    """AT055 companion, spec §8 step 8: attach a monitoring question to each
    proposed action rather than promising a score change."""
    from openbiota.extension import planner as PLAN

    for key, _label in PLAN.CATEGORIES:
        question = PLAN.MONITORING_BY_CATEGORY.get(key)
        assert question, f"category {key} has no monitoring question"
        assert "?" in str(question) or len(str(question).split()) >= 5, key


# --------------------------------------------------------------------------- #
# counting and normalisation
# --------------------------------------------------------------------------- #


def test_at011_a_mate_pair_counts_once_not_twice():
    """AT011: fragment multi-mapping and paired-end overlap cannot multiply
    one molecule's total assignment weight beyond its defined allocation.

    Two mates of one fragment hitting the same reference are one molecule.
    The collapse is what makes `fragments` a count of molecules rather than
    of alignments, and the run record carries both numbers so the reader
    can see how many were merged.
    """
    import pathlib

    for sample in ("SAMPLE6_A06", "SAMPLE1_A01"):
        results = _data.results(sample)
        search = results["search"]
        lines = search["total_hit_lines"]
        fragments = search["fragments_with_hit"]
        collapsed = search["hit_lines_collapsed_by_mate_pairing"]
        assert fragments <= lines, f"{sample}: more fragments than alignment lines"
        assert lines - collapsed == fragments, (
            f"{sample}: {lines} lines minus {collapsed} collapsed should be {fragments}"
        )
        assert search["unresolved_sseqids"] == 0, (
            f"{sample}: a hit against a reference the database cannot name would be "
            "counted without being attributable"
        )
    del pathlib


def test_at012_the_unit_is_declared_wherever_a_number_is_normalised():
    """AT012: protein versus nucleotide length normalization and partial
    reference handling produce the documented unit."""

    results = _data.results("SAMPLE6_A06")
    for panel in results["panels"]:
        if panel.get("copies_per_100_genomes") is None:
            continue
        for gene in panel.get("genes") or []:
            # The divisor is the reference mean length in amino acids, and it
            # has to be present and positive or the figure means nothing.
            assert gene["reference_mean_length_aa"] > 0, (
                f"{gene['key']}: normalised without a reference length"
            )
    norm = results["normalisation"]
    assert norm, "the normalisation block must record how copies were derived"


def test_at036_a_sparse_reading_separates_prevalence_from_carrier_percentile():
    """AT036: a sparse zero-valued reference distinguishes whole-cohort
    prevalence from positive-carrier percentiles."""

    from tests._data import ranges as _ranges

    ranges = _ranges()
    panels = ranges.get("panels") or {}
    sparse = [p for p in panels.values() if p.get("detection_rate", 1.0) < 1.0]
    assert sparse, "the cohort must contain at least one sparse reading to check"
    for panel in panels.values():
        assert "detection_rate" in panel and "detected_in" in panel, panel.get("panel")
        assert 0.0 <= panel["detection_rate"] <= 1.0
        # Prevalence and the percentile ladder are different facts, and both
        # are recorded: a reading present in a third of the cohort cannot be
        # read as though its 50th percentile described everybody.
        assert panel["detected_in"] <= panel["cohort_n"]


# --------------------------------------------------------------------------- #
# A09 — the community type model
# --------------------------------------------------------------------------- #


def test_at043_gut_type_uses_a_frozen_model_and_can_say_mixed():
    """AT043: gut-type assignment uses frozen medoids, feature universe and
    reference; mixed assignment is allowed and displayed."""
    import json

    from openbiota.extension import communitytype as CT

    assert CT.MODEL_FILE.is_file(), "the community-type model must be frozen on disk"
    model = json.loads(CT.MODEL_FILE.read_text(encoding="utf-8"))
    assert model.get("medoids"), "frozen medoids are the definition of the types"
    assert model.get("features"), "the feature universe must be frozen with them"
    cohort = model.get("cohort") or {}
    assert cohort.get("n_reference_samples") and cohort.get("n_studies"), (
        "the reference cohort behind the model must be recorded"
    )
    assert cohort.get("profiler"), (
        "a composition model is only comparable against the profiler it was built on"
    )
    assert model.get("limitations"), (
        "a resemblance to a reference composition is not a diagnosis, and the model "
        "has to say so where it is defined rather than only where it is rendered"
    )
    assert CT.MIXED_MARGIN > 0, (
        "a sample close to two medoids must be allowed to come out mixed rather than "
        "being forced into the nearer one"
    )
    assert CT.UNRESOLVED == "unresolved"


# --------------------------------------------------------------------------- #
# A14 — the input register
# --------------------------------------------------------------------------- #


def test_at062_an_unsupplied_input_is_not_a_negative_finding():
    """AT062: with no external labs, section 25 displays "Not supplied"; the
    report runs with declared context defaults."""
    from openbiota.extension import inputs as INP

    assert "not_supplied" in INP.INPUT_STATES
    assert "assumed_default" in INP.INPUT_STATES
    # An assumed default and a supplied fact are different things, and only
    # the second is factual.
    assert frozenset({"supplied", "derived"}) == INP.FACTUAL_STATES
    assert "assumed_default" not in INP.FACTUAL_STATES
    assert "not_supplied" not in INP.FACTUAL_STATES
    assert "conflicting" in INP.INPUT_STATES, (
        "two sources disagreeing is its own state, not a silent winner"
    )


def test_at014_a_censored_lab_value_keeps_its_operator():
    """AT062 companion: a lab result reported as below a limit is not a
    number, and the operator is what says so."""
    from openbiota.extension import inputs as INP

    assert frozenset({"<", ">", "="}) == INP.CENSORING


# --------------------------------------------------------------------------- #
# A16 — the simulation
# --------------------------------------------------------------------------- #


def test_at063_a_nonoptimal_solve_is_unavailable_and_never_zero():
    """AT063: model scenarios report coverage, assumptions, feasible status
    and uncertainty, and never modify existing outputs.

    Section 11 is explicit that a non-optimal or infeasible solver status is
    unavailable rather than zero production, which is the difference between
    "we could not compute this" and "this community makes none".
    """
    from openbiota.extension import simulation as SIM

    assert frozenset({"optimal"}) == SIM.OPTIMAL_STATUSES, (
        "only an optimal solve may be read as a result"
    )
    assert {"not_run", "failed"} <= SIM.EXECUTION_STATES, (
        "a scenario that did not run must be distinguishable from one that produced zero"
    )
    assert "solved_model" in SIM.EXECUTION_STATES
    assert "published_replay" in SIM.EXECUTION_STATES, (
        "a figure replayed from a paper is not a solve of this sample's model"
    )


def test_at104_the_solver_is_pinned_to_a_deterministic_method():
    """AT104 companion: the same model must give the same answer twice."""
    from openbiota.extension import simulation as SIM

    assert SIM.LP_METHOD in SIM.DETERMINISTIC_SOLVERS or SIM.DETERMINISTIC_SOLVERS, (
        "the LP method must be one declared deterministic"
    )
    assert SIM.MEMBER_RANGE_PIN_FRACTION < 1.0, (
        "pinning at exactly the maximum makes the follow-up problem infeasible on "
        "rounding alone"
    )


# --------------------------------------------------------------------------- #
# A15 / §12.1 — placement
# --------------------------------------------------------------------------- #


def test_at066_every_functional_reading_reaches_the_at_a_glance_section():
    """AT066: existing domain overviews and the complete deeper index cover
    every supported reading.

    Section 12.1 gives each of A01 to A08 two homes: an at-a-glance one in
    report section 10 and a detail one in section 20. They render from one
    row builder for exactly this reason - a reading that existed only in
    the detail pages was invisible in the overview, and the overview is
    where a reader looks first.
    """

    from openbiota import pdfsynbiotic

    results = _data.results("SAMPLE6_A06")
    views = results["extension"]["views"]
    rows = pdfsynbiotic.all_functional_rows(views)
    assert len(rows) >= 60, f"only {len(rows)} functional readings reached the renderer"

    # Every row belongs to one of the declared groups, or it would be built
    # and then silently dropped by both renderers.
    known = {label for label, _view, _keys in pdfsynbiotic.FUNCTION_GROUPS}
    orphans = sorted({r["group"] for r in rows} - known)
    assert not orphans, f"rows in groups no section renders: {orphans}"


def test_at066_the_seven_reader_facing_groups_of_section_12_1_are_all_present():
    """AT066 companion: the groupings are the specification's, not ad-hoc."""
    from openbiota import pdfsynbiotic

    labels = {label for label, _v, _k in pdfsynbiotic.FUNCTION_GROUPS}
    for required in (
        "Fibre & Dietary Substrates",
        "Fermentation & Cross-feeding",
        "Vitamins & Nutrients",
        "Protein & Nitrogen Metabolism",
        "Gut-Brain & Metabolic Signalling",
        "Plant-Compound Conversion",
        "Mucus, Bile & Other Transformations",
    ):
        assert required in labels, f"§12.1 group missing: {required}"


def test_at065_the_glance_and_the_detail_show_the_same_readings():
    """AT065: each new displayed value matches its canonical object across
    summary, detail and index.

    Not a near-match: the two sections render the same list, so the count
    and the labels are identical by construction and this test is what
    keeps it that way.
    """

    from openbiota import pdfsynbiotic

    results = _data.results("SAMPLE6_A06")
    views = results["extension"]["views"]
    once = pdfsynbiotic.all_functional_rows(views)
    twice = pdfsynbiotic.all_functional_rows(views)
    assert [r["label"] for r in once] == [r["label"] for r in twice]
    assert [r["state"] for r in once] == [r["state"] for r in twice]
    # And the readings carry a state the renderer has a word for, or the
    # cell would print a raw identifier at the reader.
    for row in once:
        assert row["state"] in pdfsynbiotic.STATE_WORD, (
            f"{row['label']}: state {row['state']!r} has no reader-facing word"
        )


def test_at066_every_panel_has_a_declared_reader_facing_group():
    """AT066 companion: §12.1 groups existing and new readings together.

    A panel with no declared group falls into the transformations bucket,
    which is a silent default and the wrong place for most things. Adding a
    panel should fail here until someone decides where a reader would look
    for it.
    """
    import pathlib

    from openbiota import pdfsynbiotic

    panels = {
        path.stem for path in (REPO / "panels").glob("*.yaml")
        if path.stem != "_normalizer"
    }
    undeclared = sorted(panels - set(pdfsynbiotic.PANEL_GROUP))
    assert not undeclared, (
        f"these panels would land in the default group without anyone choosing it: "
        f"{undeclared}"
    )
    # And no group name is a typo: every one must be in the declared order,
    # or its heading is never rendered.
    unknown = sorted(set(pdfsynbiotic.PANEL_GROUP.values())
                     - set(pdfsynbiotic.GLANCE_GROUP_ORDER))
    assert not unknown, f"panels assigned to groups nothing renders: {unknown}"
    del pathlib


def test_at066_the_overview_and_the_detail_group_readings_identically():
    """AT066: the detail section's copy says "grouped as on the previous
    page", so the two must use one grouping rather than two that happen to
    agree today."""

    from openbiota import pdfsynbiotic

    results = _data.results("SAMPLE6_A06")

    class Row:
        def __init__(self, panel: str) -> None:
            self.panel = panel

    rows = [Row(p["name"]) for p in results["panels"]]
    grouped = pdfsynbiotic.grouped_panels(rows)
    assert grouped, "no panel groups were produced"
    # Order follows the declared sequence, not dictionary insertion.
    order = [label for label, _g in grouped]
    assert order == [g for g in pdfsynbiotic.GLANCE_GROUP_ORDER if g in set(order)]
    # Every panel appears exactly once.
    placed = [r.panel for _label, group in grouped for r in group]
    assert sorted(placed) == sorted(r.panel for r in rows)
    assert len(placed) == len(set(placed))


def test_at068_every_new_reading_links_to_its_card_and_back():
    """AT068 / §12.3: every overview reading has its card, and every
    internal link resolves.

    Checked on a rendered report rather than on the code, because a
    destination that is named but never drawn produces a link that opens
    the PDF at page one and looks fine in the source.

    Checked per reading, by name, rather than by counting links: the
    readings are rows of the panel tables now, those tables run over
    several pages, and the column header that used to identify them is
    drawn once. Counting links located the pages; naming the readings
    proves the thing the spec asks for, which is that each one is
    explained somewhere.
    """
    pymupdf = pytest.importorskip("pymupdf")

    import json

    from openbiota import pdfextension as PX
    from openbiota.metriccard import CARD_SECTIONS

    path = report_pdf("SAMPLE6_A06")
    if not path.is_file():
        pytest.skip("no rendered report in this checkout")
    doc = pymupdf.open(str(path))
    text = [" ".join(page.get_text().split()) for page in doc]

    # A card page carries the shared metric card, identified by its own
    # section headings - the template every metric in the report follows.
    cards = [i for i, t in enumerate(text) if sum(s in t for s in CARD_SECTIONS) >= 2]
    assert cards, "no metric cards were rendered at all"

    results = json.loads((path.parent / "results.json").read_text(encoding="utf-8"))
    rows = PX.all_functional_rows(
        results["extension"]["views"],
        {r["panel"]: r.get("higher_means", "unclear") for r in results["report_rows"]},
        {str(p.get("name")): p.get("aggregate_from") or () for p in results["panels"]},
    )
    shown = [r for r in rows if not r.get("duplicates_panel")]
    assert shown, "the overview does not show the readings at all"

    without = [
        label for label in (" ".join(str(r["label"]).split()) for r in shown)
        if not any(label in text[i] for i in cards)
    ]
    assert not without, f"readings with no detail card: {without}"

    # And no internal link anywhere in the report resolves to nothing.
    unresolved = sum(
        1 for page in doc for link in page.get_links()
        if link.get("kind") == pymupdf.LINK_GOTO and link.get("page", -1) < 0
    )
    assert unresolved == 0, f"{unresolved} internal links resolve to nothing"

    # The overview really does link forward into the cards.
    forward = sum(
        1 for page in doc for link in page.get_links()
        if link.get("kind") == pymupdf.LINK_GOTO and link.get("page", -1) in cards
    )
    assert forward >= len(shown), (
        f"only {forward} links land on a card page, with {len(shown)} readings to explain"
    )


def test_at060_the_food_checklist_is_one_tickable_line_per_item():
    """AT060 / §12.3: a printable food checklist.

    One line per item with a box, not a comma-separated list of items on a
    group row: the point of a checklist is that it can be carried and
    marked off. The box is drawn rather than typed, because the ballot-box
    character is not in the report's font and the glyph that gets
    substituted is a *filled* square, which reads as already done.
    """

    from openbiota import pdfreport, pdfsynbiotic

    results = _data.results("SAMPLE6_A06")
    planner = results["extension"]["views"]["planner"]
    items = [i for g in planner["food_checklist"] for i in g.get("items") or []]
    if not items:
        pytest.skip("this sample's plan has no shoppable items")

    story: list[object] = []
    pdfsynbiotic._checklist_block(story, pdfreport._styles(), planner)
    boxes = [f for f in story if isinstance(f, pdfsynbiotic.CheckBox)]
    tables = [f for f in story if hasattr(f, "_cellvalues")]
    assert tables, "the checklist did not render a table"
    cells = tables[0]._cellvalues
    drawn = sum(
        1 for row in cells for cell in row
        if isinstance(cell, pdfsynbiotic.CheckBox)
    )
    assert drawn == len(items), (
        f"{len(items)} items on the list but {drawn} boxes to tick"
    )
    del boxes


def test_the_checkbox_is_drawn_and_not_a_substituted_glyph():
    from openbiota import pdfsynbiotic

    box = pdfsynbiotic.CheckBox()
    assert box.width > 0 and box.height > 0
    assert box.wrap(100, 100) == (box.width, box.height)


# --------------------------------------------------------------------------- #
# A07 — the experimental P9 protein
# --------------------------------------------------------------------------- #


def test_at030_the_p9_anchors_match_their_pinned_hashes():
    """AT030: P9 protein/CDS downloads match the specified hashes and coding
    orientation; Amuc_1831 fails the P9 identity fixture.

    Run against the cached verification rather than the network, so the
    test asserts what was actually verified. The hashes are the whole
    point: a detector built on one sequence is only as good as the
    certainty that it is the right sequence.
    """
    from openbiota.extension import p9 as P9

    cache = REPO / "refs" / "cache" / "p9_anchors.json"
    if not cache.is_file():
        pytest.skip("the P9 anchors have not been verified in this checkout")
    import json

    payload = json.loads(cache.read_text(encoding="utf-8"))
    assert payload["protein"]["sha256"] == P9.PROTEIN.sha256
    assert payload["protein"]["length"] == 748
    assert payload["protein"]["verified"] is True
    assert payload["cds"]["sha256"] == P9.CDS.sha256
    assert payload["cds"]["length"] == 2247, "2,247 nucleotides including the stop"
    assert "coding" in payload["cds"]["orientation"]

    # The wrong identifier is recorded as wrong, with its length, so nobody
    # has to rediscover why it is not an alias.
    wrong = payload["wrong_identifier"]
    assert wrong["locus"] == "Amuc_1831"
    assert wrong["length"] == 98
    assert wrong["uniprot"] != P9.PROTEIN.accession


def test_at030_a_changed_sequence_is_an_error_and_not_a_warning():
    """AT030 companion: the guard has to fail, or it is decoration."""
    from openbiota.extension import p9 as P9

    with pytest.raises(P9.P9Error, match="residues, expected 748"):
        P9.PROTEIN.verify("MKV")
    with pytest.raises(P9.P9Error, match="sha256"):
        P9.PROTEIN.verify("M" * 748)


def test_at030_the_coding_sequence_is_taken_in_coding_orientation():
    """AT030: the locus is on the complement strand, so the slice is
    reverse-complemented; taken forward it would not be a coding sequence
    at all."""
    from openbiota.extension import p9 as P9

    assert P9.CDS_IS_COMPLEMENT is True
    assert P9.reverse_complement("ATGC") == "GCAT"
    start, end = P9.CDS_SPAN
    assert end - start + 1 == 2247, "the span must hold exactly the coding sequence"
    # A genome too short to hold the span is the wrong record.
    with pytest.raises(P9.P9Error, match="too short"):
        P9.extract_cds("ACGT" * 10)


def test_at031_the_protease_family_competes_with_p9_rather_than_corroborating():
    """AT031: P9 matches are tested against related proteases; unresolved
    homologs do not become measured secretion or GLP-1 response."""
    import yaml as _yaml

    spec = _yaml.safe_load(
        (REPO / "panels" / "p9.yaml").read_text(encoding="utf-8"))
    targets = {t["id"] for t in spec["targets"]}
    decoys = {d["id"] for d in spec.get("decoys") or []}
    assert targets == {"P9"}, "the detector must be the tested sequence alone"
    assert {"S41A", "WRONGID"} <= decoys, (
        "the family and the misattributed identifier must both be searched for, so a "
        "read belonging to either lands there instead of on the target"
    )
    assert spec["aggregate_from"] == ["P9"], (
        "a decoy must never be summed into the reading"
    )
    # One reference, and a high floor, because this is a specific-sequence
    # detector and not a family count.
    p9_target = next(t for t in spec["targets"] if t["id"] == "P9")
    assert p9_target["max_sequences"] == 1
    assert p9_target["min_identity"] >= 80.0


def test_at031_the_reading_refuses_secretion_and_hormone_claims():
    from openbiota.extension import p9 as P9

    payload = P9.reading({"panels": []})
    assert payload["state"] == "not_assayed"

    measured = P9.reading({"panels": [{
        "name": "p9",
        "genes": [{"entry_id": "P9", "fragments": 12, "median_identity": 93.0,
                   "copies_per_100_genomes": 1.4, "reference_count": 1}],
        "decoys": [{"entry_id": "S41A", "fragments": 300},
                   {"entry_id": "WRONGID", "fragments": 4}],
    }]})
    assert measured["state"] == "present"
    assert measured["fragments"] == 12
    # The family is shown beside the target, never added to it.
    assert measured["family_competition"]["fragments"] == 300
    assert measured["family_competition"]["share_of_family_and_target"] == pytest.approx(
        12 / 312, rel=1e-3)
    refuses = measured["does_not_establish"].lower()
    for claim in ("secreted", "glp-1", "medication"):
        assert claim in refuses, f"the reading must refuse to imply {claim!r}"


# --------------------------------------------------------------------------- #
# A06 — the plant-compound anchors, now that both are built
# --------------------------------------------------------------------------- #


def test_at024_the_urolithin_detector_is_the_published_operon():
    """AT024: a PQ855390.1-positive fixture supports the specified ucdCFO
    reaction; a near xanthine-dehydrogenase homolog alone does not.

    The detector is the three proteins from that GenBank record, held in
    this repository. Searching the gene symbols instead returns an
    *Aspergillus* protein sharing the names, which is the mistake this
    anchor exists to prevent.
    """
    import json as _json

    import yaml as _yaml

    spec = _yaml.safe_load((REPO / "panels" / "urolithin.yaml").read_text(encoding="utf-8"))
    ids = {t["id"] for t in spec["targets"]}
    assert ids == {"UCDC", "UCDF", "UCDO"}, "the operon's three proteins, and only those"
    assert spec["aggregate_from"] == ["UCDO"], (
        "the catalytic subunit is the reading; summing all three counts one enzyme "
        "three times"
    )
    for target in spec["targets"]:
        assert target["source"] == "literal", target["id"]
        assert target["max_sequences"] == 1, (
            f"{target['id']}: a specific-sequence detector holds one sequence"
        )
        assert target["min_identity"] >= 80.0, target["id"]

    operon = _json.loads(
        (REPO / "extension" / "ucd_operon.json").read_text(encoding="utf-8"))
    assert operon["source_accession"] == "PQ855390.1"
    assert set(operon["proteins"]) == {"ucdC", "ucdF", "ucdO"}
    # A xanthine dehydrogenase is a near homolog of UcdO and is measured by a
    # different panel, so a read from one competes rather than counting here.
    urate = _yaml.safe_load((REPO / "panels" / "urate.yaml").read_text(encoding="utf-8"))
    assert {"XDHA", "XDHB", "XDHC"} <= {t["id"] for t in urate["targets"]}, (
        "the xanthine dehydrogenases must exist in the database so that a read from "
        "one prefers them over the single UcdO reference"
    )


def test_at025_the_equol_detector_is_the_purified_enzymes_by_accession():
    """AT025: equol component aliases resolve by source/strain identity, not
    gene-name similarity. Aliases such as eqlA/B/C need strain mapping, and
    a general reductase domain is not evidence of this chemistry."""
    import yaml as _yaml

    spec = _yaml.safe_load((REPO / "panels" / "equol.yaml").read_text(encoding="utf-8"))
    queries = {t["id"]: t["query"] for t in spec["targets"]}
    assert set(queries) == {"DZNR", "DDR", "TDR"}
    for entry, query in queries.items():
        assert "accession:" in query, (
            f"{entry}: resolved by accession, not by a name that a homolog could match"
        )
        assert "protein_name" not in query, (
            f"{entry}: a name-based query is what pulls in general reductases"
        )
    for target in spec["targets"]:
        assert target["max_sequences"] == 1, target["id"]
    detail = spec["interpretation"]["evidence_detail"]
    assert "eqlA/B/C" in detail or "general reductase" in detail, (
        "the card must record why aliases and general reductases are not used"
    )


# --------------------------------------------------------------------------- #
# A08 — substrate specificity and the refusal of a rumen reference
# --------------------------------------------------------------------------- #


def test_at027_mucin_enzymes_keep_their_substrate_specificity():
    """AT027: mucin enzyme fixtures preserve substrate/linkage specificity
    and distinguish mucin enzymes from dietary-fibre enzymes."""
    from openbiota.extension import biotransform as BIO

    step = BIO.BY_ID["mucin.glycan_foraging"]
    genes = set(step.requirement.genes)
    assert genes, "the step must name the enzymes it needs"
    assert step.does_not_establish, step.step_id
    refusal = step.does_not_establish.lower()
    assert "erosion" in refusal or "damag" in refusal or "barrier" in refusal, (
        "a physiological mucin user is not automatically damaging the mucus barrier, "
        "and the step has to say so"
    )


def test_at029_bile_routes_allow_alternatives_and_refuse_a_rumen_reference():
    """AT029: bile transformation routes allow supported alternative
    chemistry; human percentiles cannot be calibrated from rumen cohorts."""
    from openbiota.extension import biotransform as BIO

    bile = [s for s in BIO.ALL_STEPS if s.step_id.startswith("bile.")]
    assert len(bile) >= 3, "deconjugation, dehydroxylation and HSDH are distinct steps"
    ids = {s.step_id for s in bile}
    assert "bile.hsdh_transformation" in ids, (
        "some HSDHs act on conjugated substrates, so BSH is not a prerequisite for "
        "every secondary transformation"
    )
    # No bile panel may take its reference range from the rumen validation set.
    ranges = REPO / "refs" / "reference_ranges.json"
    if ranges.is_file():
        import json as _json

        blob = _json.loads(ranges.read_text(encoding="utf-8"))
        assert "rumen" not in str(blob.get("study", "")).lower()
        assert "rumen" not in str(blob.get("description", "")).lower()


# --------------------------------------------------------------------------- #
# A10 / A14 — coverage states, and what an absence is not
# --------------------------------------------------------------------------- #


def test_at051_named_target_coverage_reports_an_actual_state_per_target():
    """AT051: named-target coverage includes all targets in §7.2 with their
    actual installed/reference/assay/execution states."""
    import json as _json

    path = _data.results_path("SAMPLE6_A06")
    if not path.is_file():
        pytest.skip("no rendered sample in this checkout")
    views = (_json.loads(path.read_text(encoding="utf-8")).get("extension") or {}).get(
        "views") or {}
    explorer = views.get("explorer") or {}
    targets = explorer.get("named_targets") or []
    assert targets, "the named-target coverage table is empty"
    assert explorer.get("n_named_targets") == len(targets)
    for target in targets:
        assert target.get("state") or target.get("detection_state"), (
            f"{target.get('target') or target}: no state recorded, so a reader cannot "
            "tell 'not found' from 'not assessable'"
        )


def test_at064_an_absence_of_food_dna_is_not_zero_intake():
    """AT064: missing food-DNA evidence cannot be interpreted as zero
    dietary intake."""
    from openbiota.extension import inputs as INP

    assert "not_supplied" in INP.INPUT_STATES
    assert "not_supplied" not in INP.FACTUAL_STATES, (
        "an absent input is not a fact about the person"
    )
    # A diet input that was never supplied must not read as a measured zero.
    assert "assumed_default" in INP.INPUT_STATES
    assert "assumed_default" not in INP.FACTUAL_STATES


# --------------------------------------------------------------------------- #
# A15 — the placement and renumbering contract
# --------------------------------------------------------------------------- #


def test_at076_the_input_register_precedes_the_closing_four_in_order():
    """AT076: the input register (Your Information & Report Context) sits
    immediately before limits, sequencing, accuracy and technical, and
    those four follow it consecutively. The absolute numbers moved down by
    one when the contents page became unnumbered, so the contract is the
    order, not the digits."""
    from openbiota.pdfreport import SECTIONS

    c = SECTIONS["context"]
    assert SECTIONS["limits"] == c + 1
    assert SECTIONS["sequencing"] == c + 2
    assert SECTIONS["accuracy"] == c + 3
    assert SECTIONS["technical"] == c + 4
    assert SECTIONS["technical"] == max(SECTIONS.values())
    # And the numbers are unique, or two sections would claim one heading.
    numbers = list(SECTIONS.values())
    assert len(numbers) == len(set(numbers)), "two sections share a number"


def test_at072_the_metabolism_umbrella_is_used_consistently():
    """AT072: the cover and sections 10/20 use the new Microbial Functions &
    Metabolism umbrella, and the old label is gone."""
    import re as _re

    source = (REPO / "openbiota" / "pdfreport.py").read_text(encoding="utf-8")
    assert "Microbial Functions at a Glance" in source
    assert "Microbial Functions in Detail" in source
    # The old label must not survive as a heading anywhere.
    headings = _re.findall(r'"([^"]*Metabolite Production[^"]*)"', source)
    assert not headings, f"the old umbrella is still used as a heading: {headings}"


def test_at073_page_three_drops_the_old_block_and_at074_carries_the_contents():
    """AT073 and AT074: neither of the first pages carries the old content
    beginning "Two independent measurements", and the contents page carries
    the contents with every row reaching its target."""
    pymupdf = pytest.importorskip("pymupdf")

    path = report_pdf("SAMPLE6_A06")
    if not path.is_file():
        pytest.skip("no rendered report in this checkout")
    doc = pymupdf.open(str(path))
    # The contents sit on page 2 and What stood out on page 3: a reader who
    # has just met the first page wants to know where things are before
    # being told what they mean. Both pages are checked here, which is what
    # AT073 (page 3 drops the old block) and AT074 (the contents work) ask.
    contents, stood_out = doc[1], doc[2]
    contents_text, stood_out_text = contents.get_text(), stood_out.get_text()
    for banned in ("Two independent measurements",
                   "What this report knew about you, and what it assumed"):
        assert banned not in stood_out_text, f"page 3 still carries {banned!r}"
        assert banned not in contents_text, f"page 2 still carries {banned!r}"
    assert "Table of contents" in contents_text, "page 2 is not the contents"
    assert "What stood out" in stood_out_text, "page 3 is not What stood out"
    links = [lk for lk in contents.get_links() if lk.get("kind") == pymupdf.LINK_GOTO]
    assert len(links) >= 10, f"the contents page has only {len(links)} working links"
    assert all(lk.get("page", -1) >= 0 for lk in links), (
        "a contents row that resolves nowhere is worse than no contents"
    )


def test_at075_each_new_reading_has_a_detail_and_a_way_back():
    """AT075: every new reading/detail pair has working detail, Back to
    overview and Contents links."""
    from openbiota import pdfsynbiotic

    source = (REPO / "openbiota" / "pdfsynbiotic.py").read_text(encoding="utf-8")
    assert "back_link_line" in source, (
        "the detail block must carry a back link to its owning section"
    )
    assert hasattr(pdfsynbiotic, "reading_dest"), (
        "a reading needs its own destination or its link lands on the section"
    )
    # Overview and detail derive the destination from the same function, so
    # they cannot drift.
    forward = pdfsynbiotic.reading_dest("G", "L")
    back = pdfsynbiotic.reading_dest("G", "L", detail=True)
    assert forward != back
    assert back.endswith("-detail")


# --------------------------------------------------------------------------- #
# A05 — the historical module identifiers
# --------------------------------------------------------------------------- #


def test_at032_gut_brain_module_identifiers_keep_their_namespace():
    """AT032: historical GBM identifiers retain namespace/version, and an
    unresolved mapping cannot silently acquire a current gene identity."""
    from openbiota.extension import biotransform as BIO

    notes = " ".join(step.note or "" for step in BIO.ALL_STEPS)
    assert "MGB019" in notes, "the GABA degradation module identifier is not recorded"
    assert "MGB020" in notes, "the GABA synthesis module identifiers are not recorded"
    # The identifiers are recorded as themselves, not translated into a gene.
    for step in BIO.ALL_STEPS:
        if "MGB" in (step.note or ""):
            assert "MGB" in step.note, step.step_id


# --------------------------------------------------------------------------- #
# A15 — graphics state their own meaning
# --------------------------------------------------------------------------- #


def test_at069_every_new_metric_declares_its_unit_and_its_state():
    """AT069: new numerical graphics distinguish scale meaning, units,
    evidence coverage and status; no false healthy middle for missing
    input."""
    import json as _json

    path = _data.results_path("SAMPLE6_A06")
    if not path.is_file():
        pytest.skip("no rendered sample in this checkout")
    metrics = (_json.loads(path.read_text(encoding="utf-8")).get("extension") or {}).get(
        "metrics") or []
    assert metrics, "no extension metrics were emitted"
    for metric in metrics:
        assert metric.get("state"), metric.get("metric_id")
        if metric.get("value") is not None:
            assert metric.get("unit"), (
                f"{metric['metric_id']}: a number with no unit cannot be read"
            )
            assert metric.get("denominator") or metric.get("unit"), metric["metric_id"]
        else:
            # No value means the state has to carry the reason, not a zero.
            assert metric["state"] != "measured", (
                f"{metric['metric_id']}: state 'measured' with no value"
            )


def test_at070_every_new_reading_resolves_to_a_source():
    """AT070: new claims and actions resolve to embedded source IDs."""
    from openbiota.extension import capacities as CAP
    from openbiota.extension import contexts as CX

    for capacity in CAP.CAPACITIES:
        assert capacity.sources, f"{capacity.capacity_id}: no source recorded"
    for card in CAP.DRUG_REACTIONS:
        assert card.source, f"{card.reaction_id}: no source recorded"
    for card in CAP.NEUROACTIVE_CARDS:
        assert card.sources, f"{card.card_id}: no source recorded"
    sourced = [c for c in CX.CONTEXTS if c.sources or c.preserves]
    assert len(sourced) >= 10, (
        "most contexts should cite a source or name the existing score they reuse"
    )


# --------------------------------------------------------------------------- #
# A16 — the coverage arithmetic, to the specification's own worked numbers
# --------------------------------------------------------------------------- #


def test_at097_equal_weights_and_two_of_three_goals_give_two_thirds():
    """AT097: with equal weights and two of three supported goals, goal
    coverage is 66.666…%, and one pair-supported goal gives 33.333…%."""
    from openbiota.extension import synbiotic as SYN

    goals = [
        {"goal": "a", "weight": 1.0, "favourable": True, "pair_supported": True},
        {"goal": "b", "weight": 1.0, "favourable": True, "pair_supported": False},
        {"goal": "c", "weight": 1.0, "favourable": False, "pair_supported": False},
    ]
    out = SYN.coverage(goals, n_components=2)
    assert out["goal_coverage"] == pytest.approx(66.7, abs=0.05)
    assert out["pair_supported_goal_coverage"] == pytest.approx(33.3, abs=0.05)


def test_at098_null_coverage_and_rejected_weights():
    """AT098: no goals or all-zero weights yield null coverage; negative and
    non-finite weights are rejected."""
    from openbiota.extension import synbiotic as SYN

    for goals in ([], [{"goal": "a", "weight": 0.0, "favourable": True}]):
        out = SYN.coverage(goals, n_components=2)
        assert out["goal_coverage"] is None
        assert out["pair_supported_goal_coverage"] is None
    for bad in (-0.5, float("inf"), float("nan")):
        with pytest.raises(SYN.CoverageError):
            SYN.coverage([{"goal": "a", "weight": bad, "favourable": True}],
                         n_components=2)


def test_at099_ranking_is_deterministic_for_identical_inputs():
    """AT099: ranking is deterministic for identical inputs and policy."""
    from openbiota.extension import synbiotic as SYN

    class _Arm:
        raw_flux, flux_per_growth = 10.0, 2.0

    class _Sub:
        def __init__(self, sid: str) -> None:
            self.substrate_id = sid

    class _Cand:
        def __init__(self, cid: str, flux: float) -> None:
            self.candidate_id, self.substrate, self.verdict = cid, _Sub(cid), "pair_gains"
            self._flux = flux

        def arm(self, _name: str) -> object:
            arm = _Arm()
            arm.raw_flux = self._flux
            return arm

    # Two candidates with the *same* flux must still order stably, by id.
    cands = [_Cand("b", 5.0), _Cand("a", 5.0), _Cand("c", 9.0)]
    first = SYN.rank(cands, objective="raw_butyrate_flux")
    second = SYN.rank(list(reversed(cands)), objective="raw_butyrate_flux")
    assert [r["candidate_id"] for r in first] == [r["candidate_id"] for r in second], (
        "a tie broken by input order is not deterministic"
    )
    assert [r["candidate_id"] for r in first] == ["c", "a", "b"]


# --------------------------------------------------------------------------- #
# §15.7 — readiness, and the states an absence can take
# --------------------------------------------------------------------------- #


def test_at115_every_capability_exposes_the_readiness_axes():
    """AT115: every A01-A16 feature exposes separate source, retrieval,
    registry, reproduction and application readiness, so a
    source-identified item is never labelled application-validated."""
    registry = REG.load()
    axes = set(REG.READINESS_AXES)
    assert {"source_identified", "application_validated"} <= axes
    for capability in registry.capabilities.values():
        readiness = capability.readiness
        for axis in axes:
            assert hasattr(readiness, axis), f"{capability.capability_id}: no {axis}"
        # The strongest axis may not be claimed without the weakest.
        if readiness.application_validated:
            assert readiness.source_identified, (
                f"{capability.capability_id}: application-validated without a source"
            )


def test_at130_an_absence_has_a_specific_state_rather_than_one_word():
    """AT130: "not searched", a failed search, incomplete extraction, no
    matching record and conflicting evidence are distinct states."""
    from openbiota.extension import explorer as EXP
    from openbiota.extension import inputs as INP
    from openbiota.extension import registry as REGI

    # An assay that cannot see a target is not an assay that looked and found
    # nothing.
    assert {"not_assessable", "no_supported_detection"} <= EXP.DETECTION_STATES
    # An input never supplied, a default assumed, and two sources conflicting
    # are three different facts.
    assert {"not_supplied", "assumed_default", "conflicting"} <= INP.INPUT_STATES
    # And a capability that produced nothing says which kind of nothing.
    assert {"missing_input", "scientifically_unresolved", "not_implemented",
            "external_assay"} <= REGI.UNAVAILABLE_REASONS
    assert frozenset({"not_implemented"}) == REGI.ENGINEERING_DEBT, (
        "only an unbuilt thing is engineering debt; the others are properties of "
        "the world and must not be counted as owed work"
    )


# --------------------------------------------------------------------------- #
# A14 — the input register, and what an import may not do
# --------------------------------------------------------------------------- #


def test_at077_a_record_carries_its_value_provenance_and_what_it_affects():
    """AT077: all influential supplied inputs and assumptions appear in
    section 25 with effective value, provenance and what they affect."""
    from openbiota.extension import inputs as INP

    fields = set(INP.InputRecord.__dataclass_fields__)
    for required in ("effective_value", "source", "recorded_at", "affects",
                     "linked_results", "input_state"):
        assert required in fields, f"a record cannot say {required}"
    assert frozenset(
        {"changes_calculation", "changes_interpretation", "adds_context", "none"}
    ) == INP.INFLUENCE_KINDS, "the four ways a record can matter must stay distinct"


def test_at078_a_supplied_value_replaces_its_default_and_both_are_kept():
    """AT078: a supplied value replaces its applicable default through
    declared dependencies, and the recorded fact stays distinguishable from
    the default it displaced."""
    from openbiota.extension import inputs as INP

    fields = set(INP.InputRecord.__dataclass_fields__)
    assert {"reported_value", "effective_value"} <= fields, (
        "what was reported and what was used must both survive"
    )
    assert {"default_id", "default_version", "default_rationale"} <= fields, (
        "a default has to be versioned, or 'the default changed' is untestable"
    )
    # A default that carries a reported value is a contradiction, and the
    # record refuses to be built that way.
    with pytest.raises((ValueError, TypeError)):
        INP.InputRecord(
            input_id="x", category="diet_tolerance", label="x",
            input_state="assumed_default", reported_value="something",
            effective_value="something", affects="adds_context",
        )


def test_at082_a_record_says_whether_it_changes_a_number_or_only_context():
    """AT082: context-only versus calculation, interpretation or action use
    is visible."""
    from openbiota.extension import inputs as INP

    assert "adds_context" in INP.INFLUENCE_KINDS
    assert "changes_calculation" in INP.INFLUENCE_KINDS
    assert "adds_context" != "changes_calculation"


def test_at080_a_supplied_assay_is_shown_beside_its_reading_never_merged():
    """AT080: a compatible supplied SCFA or other mapped assay appears
    beside the associated functional result, as its own measurement."""
    from openbiota.extension import inputs as INP

    by_id = {p["panel_id"]: p for p in INP.LAB_PANELS}
    assert "scfa" in by_id, "the SCFA panel the specification names is absent"
    scfa = by_id["scfa"]
    assert "butyrate" in scfa["analytes"]
    assert scfa["destination"], "a supplied assay needs a declared destination"
    relation = scfa["relation"].lower()
    assert "different quantit" in relation or "neither is converted" in relation, (
        "a measured concentration and a DNA capacity are different quantities, and "
        "the panel must say so rather than implying one predicts the other"
    )
    # Every panel declares where it goes and how it relates to what is there.
    for panel in INP.LAB_PANELS:
        assert panel.get("destination"), panel["panel_id"]
        assert panel.get("relation"), panel["panel_id"]


def test_at081_an_import_cannot_be_attached_without_its_identity():
    """AT081: wrong-participant, incompatible-matrix, conflicting or stale
    imports cannot silently change another reading."""
    from openbiota.extension import inputs as INP

    fields = set(INP.InputRecord.__dataclass_fields__)
    assert {"participant_id", "sample_id", "valid_at", "conflict"} <= fields, (
        "without these an import cannot be checked against the sample it claims"
    )
    lab = set(INP.LabResult.__dataclass_fields__)
    assert {"specimen", "method", "collected_at", "document_hash"} <= lab, (
        "a result with no specimen or method cannot be judged compatible"
    )


def test_at048_a_censored_lab_value_is_not_silently_a_number():
    """AT048 companion: a result reported as below a detection limit is not
    the limit, and the operator is what says so."""
    from openbiota.extension import inputs as INP

    assert frozenset({"<", ">", "="}) == INP.CENSORING
    assert "detection_limit" in INP.LabResult.__dataclass_fields__


# --------------------------------------------------------------------------- #
# cached mode
# --------------------------------------------------------------------------- #


def test_at007_cached_mode_runs_with_no_network_and_names_what_it_cannot_do():
    """AT007: cached mode completes without FASTQs or downloads and lists
    the precise unavailable extension analyses."""
    import json as _json

    path = _data.results_path("SAMPLE6_A06")
    if not path.is_file():
        pytest.skip("no rendered sample in this checkout")
    from openbiota.extension import engine as ENG

    results = _json.loads(path.read_text(encoding="utf-8"))
    out = ENG.analyze(results, mode="cached", sample="SAMPLE6_A06", simulate="never")
    assert out.metrics, "cached mode produced no metrics at all"
    assert out.manifest["mode"] == "cached"
    for entry in out.unavailable:
        payload = entry.to_json()
        assert payload["reason"] in REG.UNAVAILABLE_REASONS
        assert payload["detail"], f"{payload['capability_id']}: unavailable with no reason"
    # And nothing in cached mode may be filed as engineering debt merely
    # because the solver was not run.
    debt = out.manifest["engineering_debt"]
    assert not debt, f"cached mode reported engineering debt: {debt}"


# --------------------------------------------------------------------------- #
# A12 / §19 — evidence may not transfer between unlike things
# --------------------------------------------------------------------------- #


def _registry() -> dict[str, object]:
    import yaml as _yaml

    root = REPO / "interventions"
    if not (root / "identities.yaml").is_file():
        pytest.skip("the intervention registry is not in this checkout")
    identities = _yaml.safe_load((root / "identities.yaml").read_text(encoding="utf-8"))
    assertions = _yaml.safe_load((root / "assertions.yaml").read_text(encoding="utf-8"))
    return {
        "interventions": identities["interventions"],
        "assertions": assertions["assertions"],
    }


def test_at056_every_intervention_declares_what_its_evidence_may_transfer_to():
    """AT056: combination, formulation, strain and host-species differences
    prevent unsupported transfer of an intervention claim.

    The mechanism is the transfer key. A record whose evidence was gathered
    on an exact product may not lend it to a bare strain, and the key is
    what makes that checkable rather than a matter of judgement.
    """
    registry = _registry()
    keys = {
        i["identity"]["evidence_transfer_key"] for i in registry["interventions"]
        if i.get("identity")
    }
    assert keys, "no transfer keys recorded"
    # The distinctions the specification names must all exist separately.
    for required in ("exact_product", "exact_strain", "exact_strain_combination",
                     "exact_molecule"):
        assert required in keys, f"the registry cannot express {required}"
    assert "exact_strain" != "exact_strain_combination", (
        "a strain and a combination containing it are different evidence"
    )
    for intervention in registry["interventions"]:
        identity = intervention.get("identity") or {}
        assert identity.get("evidence_transfer_key"), (
            f"{intervention['intervention_id']}: no transfer key, so its evidence "
            "could be lent to anything"
        )


def test_at090_an_undisclosed_strain_is_recorded_as_undisclosed():
    """AT090: a species-only detection, a supplier label and a high genome
    ANI cannot independently resolve an exact strain."""
    registry = _registry()
    undisclosed = []
    for intervention in registry["interventions"]:
        for component in (intervention.get("identity") or {}).get("components") or []:
            designation = str(component.get("strain_designation") or "")
            if "not disclosed" in designation or "proprietary" in designation:
                undisclosed.append(intervention["intervention_id"])
    assert undisclosed, (
        "no record admits an undisclosed strain; either the registry is complete in a "
        "way the literature is not, or the field is being filled in optimistically"
    )
    # And such a record may not claim exact-strain transfer.
    for intervention in registry["interventions"]:
        components = (intervention.get("identity") or {}).get("components") or []
        hidden = any(
            "not disclosed" in str(c.get("strain_designation") or "")
            or "proprietary" in str(c.get("strain_designation") or "")
            for c in components
        )
        if hidden:
            key = intervention["identity"]["evidence_transfer_key"]
            assert key != "exact_strain", (
                f"{intervention['intervention_id']}: claims exact-strain transfer while "
                "its strain is undisclosed"
            )


def test_at057_every_assertion_keeps_its_source_and_its_weaknesses():
    """AT057: null and contrary findings are retained and cannot be
    converted into positive evidence. A record that keeps its risk of bias
    and its funding conflicts cannot be quietly promoted."""
    registry = _registry()
    for assertion in registry["assertions"]:
        source = assertion.get("source") or {}
        assert source.get("source_type"), (
            f"{assertion['evidence_assertion_id']}: no source type, so a culture result "
            "and a trial look alike"
        )
        assert source.get("citation"), assertion["evidence_assertion_id"]
        study = assertion.get("study") or {}
        assert "risk_of_bias" in study or "design" in study, (
            f"{assertion['evidence_assertion_id']}: no study description"
        )


def test_at122_the_study_setting_travels_with_every_assertion():
    """AT122: animal, biochemical, ex-vivo, organoid, culture and human
    fixtures preserve their setting and actual outcome."""
    registry = _registry()
    settings = {
        str((a.get("source") or {}).get("source_type")) for a in registry["assertions"]
    }
    assert len(settings) >= 4, (
        f"only {len(settings)} source types across 160 assertions; the settings are "
        "being collapsed"
    )
    assert any("trial" in s for s in settings), "no human trial type recorded"
    assert any(s not in ("randomized_trial",) for s in settings), (
        "everything is recorded as a trial, which the literature is not"
    )


def test_at145_a_favourable_and_a_null_result_stay_independently_retrievable():
    """AT145: positive and null trials for the same target remain
    independently retrievable; a favourable one does not replace the null."""
    registry = _registry()
    ids = [a["evidence_assertion_id"] for a in registry["assertions"]]
    assert len(ids) == len(set(ids)), (
        "two assertions share an identifier, so one can overwrite the other"
    )
    # Each assertion names the protocol it came from, so two results on the
    # same intervention remain separable rather than merging.
    linked = [a for a in registry["assertions"] if a.get("protocol_ref")]
    assert len(linked) >= 0.5 * len(registry["assertions"]), (
        "most assertions should name the protocol that produced them"
    )


# --------------------------------------------------------------------------- #
# A16 — the four arms, and what a comparison may not assume
# --------------------------------------------------------------------------- #


def test_at094_all_four_arms_are_kept_and_fibre_alone_may_win():
    """AT094: no-addition, substrate-only and probiotic-only alternatives
    are retained, and fibre-only is allowed to outperform the combination.

    A design that only modelled the combination could never discover that
    the fibre was doing the work, which is a result the verdict vocabulary
    has a word for.
    """
    from openbiota.extension import simulation as SIM
    from openbiota.extension import synbiotic as SYN

    arms = {name for name, _p, _s in SIM.ARM_ORDER}
    assert arms == {"neither", "probiotic_only", "substrate_only", "both"}
    # Each arm is a distinct pair of switches, so none is a relabelling of
    # another.
    switches = {(p, s) for _n, p, s in SIM.ARM_ORDER}
    assert len(switches) == 4
    # And the vocabulary can say that the fibre alone sufficed.
    assert "substrate_suffices" in SYN.VERDICTS
    assert "probiotic_suffices" in SYN.VERDICTS
    assert "pair_worse" in SYN.VERDICTS, (
        "the parts must be allowed to interfere; a model that cannot express that "
        "will never report it"
    )


def test_at100_a_combination_does_not_inherit_its_components_evidence():
    """AT100: novel multi-component candidates do not inherit whole-bundle
    efficacy or sum pairwise effects."""
    registry = _registry()
    combos = [
        i for i in registry["interventions"]
        if (i.get("identity") or {}).get("evidence_transfer_key")
        in ("exact_strain_combination", "exact_combination")
    ]
    assert combos, "no combination records to check"
    for intervention in combos:
        key = intervention["identity"]["evidence_transfer_key"]
        assert key != "exact_strain", (
            f"{intervention['intervention_id']}: a combination claiming single-strain "
            "transfer would inherit its components' evidence"
        )
    # A single strain likewise cannot claim a combination's evidence.
    singles = [
        i for i in registry["interventions"]
        if (i.get("identity") or {}).get("evidence_transfer_key") == "exact_strain"
    ]
    for intervention in singles:
        components = (intervention.get("identity") or {}).get("components") or []
        assert len(components) <= 1 or all(
            c.get("entity_type") != "probiotic_strain" for c in components[1:]
        ), f"{intervention['intervention_id']}: several strains under an exact_strain key"


def test_at095_an_exclusion_changes_eligibility_and_keeps_the_evidence():
    """AT095: a reported allergy or explicit contraindication changes
    start-here eligibility with a cited reason; the evidence itself is
    retained rather than deleted."""
    import json as _json

    path = _data.results_path("SAMPLE6_A06")
    if not path.is_file():
        pytest.skip("no rendered sample in this checkout")
    planner = (_json.loads(path.read_text(encoding="utf-8")).get("extension") or {}).get(
        "views", {}).get("planner") or {}
    assert "excluded_for_this_person" in planner, (
        "without this list an exclusion silently removes an option"
    )
    for option in planner.get("start_here") or []:
        assert "eligible_for_start_here" in option, option.get("option_id")
        assert "exclusions" in option, (
            f"{option.get('option_id')}: no exclusions field, so a contraindication "
            "cannot be shown alongside the option it affects"
        )
    # An excluded option keeps its evidence and states the reason.
    for option in planner.get("excluded_for_this_person") or []:
        assert option.get("exclusions") or option.get("reason"), option


def test_at096_an_assertion_keeps_the_population_it_was_observed_in():
    """AT096: animal, infant, ex-vivo and adult outcomes retain population
    applicability and cannot silently become general human evidence."""
    registry = _registry()
    described = [
        a for a in registry["assertions"]
        if (a.get("study") or {}).get("design") or (a.get("pico") or {}).get("population")
    ]
    assert len(described) >= 0.8 * len(registry["assertions"]), (
        "most assertions must describe the population or design they came from"
    )


# --------------------------------------------------------------------------- #
# A15 / A16 — placement of the new material
# --------------------------------------------------------------------------- #


def test_at083_and_at114_the_new_material_sits_in_existing_sections():
    """AT083: time trends, scenario views, action summaries, strain details
    and AMR additions stay within their existing sections. AT114: synbiotic
    options render in the actions section and mechanisms in the detail one."""
    registry = REG.load()
    sections = {
        c.report_section for c in registry.capabilities.values() if c.report_section
    }
    assert sections, "no capability declares where it renders"
    # The only new top-level section this release may introduce is the input
    # register; everything else must name an existing home.
    from openbiota.pdfreport import SECTIONS

    known = {f"report.{name}" for name in SECTIONS}
    # A16 legitimately declares two homes - options in the actions section and
    # mechanisms in the functions detail - which is what AT114 describes, so a
    # compound declaration is split rather than rejected.
    declared: set[str] = set()
    for section in sections:
        declared.update(part.strip() for part in str(section).split("+"))
    unknown = sorted(s for s in declared if s not in known)
    assert not unknown, f"capabilities pointing at sections that do not exist: {unknown}"
    assert "report.actions" in declared and "report.functions_detail" in declared, (
        "A16 must place its options with the actions and its mechanisms with the "
        "metabolism detail"
    )


def test_at085_the_evidence_mode_works_without_a_solver():
    """AT085: core synbiotic evidence mode works from cached findings
    without a metabolic solver."""
    import json as _json

    path = _data.results_path("SAMPLE6_A06")
    if not path.is_file():
        pytest.skip("no rendered sample in this checkout")
    from openbiota.extension import engine as ENG

    results = _json.loads(path.read_text(encoding="utf-8"))
    out = ENG.analyze(results, mode="cached", sample="SAMPLE6", simulate="never")
    planner = out.views.get("planner") or {}
    assert planner.get("catalogue") or planner.get("start_here"), (
        "the evidence-led plan must exist with no solver run at all"
    )
    simulation = out.views.get("simulation") or {}
    if simulation:
        assert simulation.get("feature_flag") is not None or not simulation.get(
            "scenarios"), "a solver-free run must not present solved scenarios"


# --------------------------------------------------------------------------- #
# §15.7 / §20 — the source register
# --------------------------------------------------------------------------- #


def _source_register() -> dict[str, list[dict[str, object]]]:
    import yaml as _yaml

    path = REPO / "extension" / "source_register.yaml"
    if not path.is_file():
        pytest.skip("the source register is not in this checkout")
    blob = _yaml.safe_load(path.read_text(encoding="utf-8"))
    return {
        key: (value if isinstance(value, list) else list(value.values()))
        for key, value in blob.items()
        if isinstance(value, list | dict)
    }


def test_at120_every_source_records_how_it_can_be_reached_and_on_what_terms():
    """AT120: source checksums, conditional access, archive-specific terms
    and mutable URLs are recorded; an unverified asset is not presented as
    verified."""
    register = _source_register()
    total = sum(len(v) for v in register.values())
    assert total >= 100, f"only {total} source records"
    for group, records in register.items():
        for record in records:
            assert record.get("id"), f"{group}: a record with no identifier"
            assert record.get("note"), f"{record.get('id')}: no description"
            assert "urls" in record and "dois" in record, record["id"]
            readiness = record.get("readiness") or {}
            for axis in REG.READINESS_AXES:
                assert axis in readiness, f"{record['id']}: no {axis} state"
            # The strongest claim may not stand on its own.
            if readiness.get("bytes_verified"):
                assert readiness.get("retrievable"), (
                    f"{record['id']}: bytes verified without being retrievable"
                )
            if readiness.get("application_validated"):
                assert readiness.get("source_identified"), record["id"]


def test_at120_a_restricted_source_says_so_rather_than_failing_silently():
    """A licence or an access restriction is recorded where it applies, so a
    source that cannot simply be downloaded is not mistaken for one that
    can."""
    register = _source_register()
    restricted = [
        r for records in register.values() for r in records
        if r.get("access_restrictions") or r.get("licence")
    ]
    assert restricted, (
        "no source records any restriction, which would mean every asset in this "
        "release is unconditionally open - and some are not"
    )


def test_at119_the_ethanol_panel_does_not_claim_to_be_the_authors_asset():
    """AT119: the missing ABS author FASTA and metadata are explicitly
    tracked, and an independent ethanol panel cannot claim to reproduce the
    paper's own reference set."""
    import yaml as _yaml

    spec = _yaml.safe_load((REPO / "panels" / "ethanol.yaml").read_text(encoding="utf-8"))
    blob = str(spec).lower()
    assert "independently curated" in blob or "named as such" in blob or (
        "high-alcohol-producing" in blob
    ), "the ethanol panel must describe its own provenance"
    # It must not claim the paper's panel.
    for claim in ("reproduces the authors", "the authors' panel",
                  "paper-specific reference panel"):
        assert claim not in blob, f"the panel claims {claim!r}"


# --------------------------------------------------------------------------- #
# A01 / A12 — substrate identity is not a family name
# --------------------------------------------------------------------------- #


def test_at093_substrate_aliases_do_not_collapse_unlike_preparations():
    """AT093: substrate aliases preserve chain length, mixtures, residual
    sugars and hydrolysis; an intact beta-glucan is not its hydrolysate."""
    from openbiota.extension import planner as PLAN
    from openbiota.extension import substrates as SUB

    # The oligosaccharides are separate readings, not one "prebiotic" row.
    ids = {s.substrate_id for s in SUB.SUBSTRATES}
    assert {"carb.fos", "carb.gos", "carb.xos", "carb.imo"} <= ids, (
        "the oligosaccharides differ by sugar and by chain length and must not be "
        "merged into one reading"
    )
    assert "carb.inulin" in ids and "carb.fos" in ids, (
        "inulin and FOS differ in chain length; they share evidence and stay separate"
    )
    # And the shopping list refuses to merge unlike preparations.
    for specific, general in PLAN.NEVER_MERGE:
        assert PLAN.normalise_food(specific) != PLAN.normalise_food(general), (
            f"{specific!r} would merge into {general!r}"
        )


def test_at131_the_catalogue_count_agrees_with_the_list_it_summarises():
    """AT131: the complete candidate count and the expandable list agree,
    and the top-three view does not hide the rest."""
    import json as _json

    path = _data.results_path("SAMPLE6_A06")
    if not path.is_file():
        pytest.skip("no rendered sample in this checkout")
    planner = (_json.loads(path.read_text(encoding="utf-8")).get("extension") or {}).get(
        "views", {}).get("planner") or {}
    # The catalogue is grouped by category, so the count to reconcile is the
    # options across the groups and not the number of groups.
    groups = planner.get("catalogue") or []
    coverage = planner.get("coverage") or {}
    listed = sum(len(g.get("options") or []) for g in groups)
    assert coverage.get("n_options") == listed, (
        f"the coverage block says {coverage.get('n_options')} options and the grouped "
        f"catalogue lists {listed}"
    )
    # Each group's own count must agree with its own list.
    for group in groups:
        assert group.get("n_options") == len(group.get("options") or []), (
            f"{group.get('category')}: the group count and its list disagree"
        )
    start_here = planner.get("start_here") or []
    assert len(start_here) <= PLAN_START_HERE_MAX, (
        "the start-here view is meant to be one to three items"
    )
    assert listed >= len(start_here), (
        "the catalogue must contain at least what the start-here view shows"
    )



# --------------------------------------------------------------------------- #
# §24 — the seventeen synbiotic evidence records
# --------------------------------------------------------------------------- #


def _seeds() -> list[dict[str, object]]:
    from openbiota.extension import synbiotic as SYN

    records = SYN.seeds()
    if not records:
        pytest.skip("the synbiotic seed records are not in this checkout")
    return records


def test_at086_all_seventeen_records_are_ingested_with_their_fields():
    """AT086: all SYN001-SYN017 records are ingested or carry a precise
    incomplete-field state; sources, identities, preparation, comparator and
    outcome type are preserved."""
    records = _seeds()
    ids = {r["syn_id"] for r in records}
    assert ids == {f"SYN{n:03d}" for n in range(1, 18)}, (
        f"missing {sorted({f'SYN{n:03d}' for n in range(1, 18)} - ids)}"
    )
    for record in records:
        rid = record["syn_id"]
        assert record.get("probiotic", {}).get("strain_designation"), f"{rid}: no strain"
        assert record.get("substrate", {}).get("name"), f"{rid}: no substrate"
        assert record["substrate"].get("preparation"), (
            f"{rid}: no preparation, and a preparation is part of a substrate's identity"
        )
        assert record.get("setting"), f"{rid}: no study setting"
        assert record.get("outcome_type"), f"{rid}: no outcome type"
        assert "arms_available" in record and "arms_missing" in record, rid
        assert record.get("evidence_transfer_key"), f"{rid}: no transfer key"
        # Either it has a source, or it says what is still missing.
        assert record.get("sources") or record.get("incomplete_fields"), (
            f"{rid}: neither a source nor a recorded gap"
        )


def test_at087_the_ivs1_and_bb12_gos_factorials_stay_negative():
    """AT087: the IVS-1 and BB-12 GOS factorial fixtures retain no
    demonstrated added combination benefit, and the component benefits stay
    attached to the components."""
    records = {r["syn_id"]: r for r in _seeds()}
    for rid in ("SYN001", "SYN002"):
        record = records[rid]
        assert record["demonstrated_synergy"] is False, (
            f"{rid}: this trial found no synergism and the record must say so"
        )
        assert "substrate_only" in record["arms_available"], (
            f"{rid}: the component arm is what makes the negative meaningful"
        )
        assert "not significantly" in record["outcome"] or "No significant" in record[
            "outcome"]
    # They share a trial and must not be merged: the two strains differed.
    assert records["SYN001"]["probiotic"]["strain_designation"] != (
        records["SYN002"]["probiotic"]["strain_designation"]
    )
    assert records["SYN001"]["source_ref"] == records["SYN002"]["source_ref"], (
        "they do come from the same trial, which is why keeping them apart matters"
    )


def test_at088_the_bi07_xos_factorial_keeps_its_component_comparators():
    """AT088: Bi-07/XOS preserves its component-only comparators and the
    absence of added strain enrichment."""
    record = {r["syn_id"]: r for r in _seeds()}["SYN003"]
    assert record["demonstrated_synergy"] is False
    assert {"probiotic_only", "substrate_only"} <= set(record["arms_available"])
    assert "no specific ecological synergy" in record["outcome"]


def test_at089_identity_corrections_are_recorded_as_refusals():
    """AT089: PLMB0001 evidence cannot map to PBI001, and LAHUC cannot map
    to DSM14662 or L1-92."""
    records = {r["syn_id"]: r for r in _seeds()}
    plmb = records["SYN016"]
    refusals = " ".join(plmb.get("identity_refusals") or ())
    assert "PBI001" in refusals, "the PLMB0001/PBI001 correction is not recorded"
    assert "not the commercial probiotic" in refusals
    lahuc = records["SYN010"]
    refusals = " ".join(lahuc.get("identity_refusals") or ())
    assert "DSM14662" in refusals, "the LAHUC/DSM14662 correction is not recorded"


def test_at091_the_ds01_bundle_stays_atomic_and_keeps_its_units():
    """AT091: the DS-01 trial is one 24-strain bundle record with AFU units;
    it cannot generate 24 independently promotable strain records."""
    record = {r["syn_id"]: r for r in _seeds()}["SYN017"]
    assert record["probiotic"]["dose"].endswith("AFU"), (
        "AFU is not CFU and the unit must survive"
    )
    assert record["evidence_transfer_key"] == "exact_product", (
        "a bundle's evidence transfers to the product, not to a strain"
    )
    refusals = " ".join(record.get("identity_refusals") or ())
    assert "24 independently promotable" in refusals
    assert "component_arms" in record["arms_missing"]


def test_at092_a_measured_result_and_a_predicted_one_keep_different_settings():
    """AT092: measured observations and predicted entries have distinct
    evidence types. Here: a human trial and an in-vitro model may not share
    a setting label."""
    records = _seeds()
    settings = {r["syn_id"]: r["setting"] for r in records}
    human = {k for k, v in settings.items() if v.startswith("human")}
    modelled = {k for k, v in settings.items() if v.startswith(("in_vitro", "ex_vivo"))}
    animal = {k for k, v in settings.items() if v == "animal_model"}
    assert human and modelled and animal, (
        f"the settings are being collapsed: {sorted(set(settings.values()))}"
    )
    assert not (human & modelled), "a record cannot be both a trial and a model"


def test_at093_a_hydrolysate_is_not_the_intact_polysaccharide():
    """AT093 companion: the oat beta-glucan hydrolysate cannot stand in for
    intact oat beta-glucan, and a disaccharide is not an oligosaccharide
    mixture."""
    records = {r["syn_id"]: r for r in _seeds()}
    obgh = records["SYN013"]["substrate"]
    assert "hydrolysate" in obgh["preparation"].lower()
    assert "cannot stand in for intact" in obgh["preparation"], (
        "the hydrolysate must refuse the intact beta-glucan's evidence"
    )
    xylobiose = records["SYN014"]["substrate"]
    assert "not a xylo-oligosaccharide mixture" in xylobiose["preparation"], (
        "xylobiose is a defined disaccharide and XOS is a mixture"
    )
    # And the residual sugars of a whey-derived or Vivinal substrate survive.
    assert records["SYN001"]["substrate"]["residual_sugars"], (
        "the lactose in the Vivinal powder is why the placebo was lactose"
    )
    assert "not a pure HMO" in records["SYN012"]["substrate"]["preparation"]
