"""Command line interface. Orchestration only — no science lives here."""

from __future__ import annotations

import argparse
import concurrent.futures
import dataclasses
import json
import os
import shutil
import subprocess
import sys
import textwrap
import time
import traceback
from collections.abc import Sequence
from pathlib import Path
from typing import Any, Final

from openbiota import __version__
from openbiota.actions import ActionPlan, build_action_plan
from openbiota.age import AgeResult, adjust_for_community_health, predict_age
from openbiota.age import bundle_dir as age_bundle_dir
from openbiota.age import load_bundle as load_age_bundle
from openbiota.batch import (
    SampleOutcome,
    outcome_from_results,
    print_comparison,
    write_comparison,
)
from openbiota.cohort import (
    DEFAULT_SEED,
    DEFAULT_STUDY,
    SAMPLING_MODES,
    TRANSPORTS,
    CohortRun,
    ReferenceRanges,
    build_ranges,
    download_cohort,
    list_runs,
    study_title,
)
from openbiota.context import resolve_context
from openbiota.depth import (
    DEFAULT_FRACTIONS,
    DepthPoint,
    DepthReport,
    assess_convergence,
    print_depth_summary,
)
from openbiota.engines.metaphlan import locate_tools
from openbiota.engines.metaphlan4 import (
    ExtendedCatalogue,
    locate_tools4,
    run_metaphlan4,
    run_sample2markers,
)
from openbiota.errors import DependencyError, OpenBiotaError
from openbiota.interpret import build_rows, cohort_note, load_optional_json, reference_json
from openbiota.logging_util import Reporter, Style, human_bytes, human_duration
from openbiota.panels import Panel, PanelSet, load_panel_set
from openbiota.pathogens.catalog import load_catalog as load_pathogen_catalog
from openbiota.pathogens.detect import PathogenBranchResult, run_pathogen_branch
from openbiota.pdfreport import build_pdf
from openbiota.preprocess import (
    MIN_USABLE_NONHOST_PAIRS,
    FastpResult,
    HeaderProvenance,
    assemble_gates,
    locate_fastp,
    read_header_provenance,
    run_fastp,
)
from openbiota.profiles import load_profile_set
from openbiota.profilevalidate import (
    MIN_PER_ARM,
    SpecificityMatrix,
    evaluate_auc,
    score_labelled_cohort,
    spread_verdict,
)
from openbiota.qc import DEFAULT_SAMPLE_READS, SampleQC, profile_sample_cached
from openbiota.refcohort import (
    DEFAULT_SNAPSHOT,
    ReferenceCohort,
    build_cohort,
    build_taxon_references,
    stability_curve,
)
from openbiota.references import (
    OUT_OF_SCOPE_SUFFIX,
    build_reference_database,
    database_fingerprint,
    diamond_version,
)
from openbiota.report import (
    build_results_json,
    render_compact,
    render_summary,
    run_diagnostics,
    write_outputs,
    write_run_log,
)
from openbiota.search import (
    MateInput,
    SampleInput,
    SearchConfig,
    check_platform,
    discover_sample,
    discover_samples,
    output_fields,
    platform_info,
    run_search,
)
from openbiota.similarity import (
    SimilarityStage,
    Subject,
    ensure_host_index,
    run_similarity_stage,
    run_taxonomic_engine,
)
from openbiota.tally import iter_lines, tally, write_panel_hit_views
from openbiota.taxongroups import load_taxon_group_set
from openbiota.taxonomy import build_profile

REPO_ROOT: Final = Path(__file__).resolve().parent.parent
SUBCOMMANDS: Final = (
    "run", "panels", "build-db", "doctor", "cohort", "validate", "depth-check",
    "taxonomic-cohort", "profiles", "validate-profiles", "build-host-index", "run-all",
    "prune", "age-model", "fmt", "biofilm", "mycobiome",
)

SENSITIVITIES: Final = (
    "default", "faster", "fast", "mid-sensitive", "sensitive",
    "more-sensitive", "very-sensitive", "ultra-sensitive",
)


# --------------------------------------------------------------------------- #
# argument parsing
# --------------------------------------------------------------------------- #


def _default_threads() -> int:
    return os.cpu_count() or 4


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="openbiota",
        description=(
            "OpenBiota Gut Health Test. Turns raw shotgun stool metagenome reads "
            "(FASTQ) into a gut-health report: metabolite pathway gene capacity normalised "
            "to bacterial genome equivalents via rpoB, species-level community "
            "composition, a published gut-health index, and resemblance to published "
            "disease-associated microbiome patterns. Research use; not a diagnostic test."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "examples:\n"
            "  openbiota                                 # full run on ./fastq, all panels\n"
            "  openbiota run --subsample 100000          # fast smoke test\n"
            "  openbiota run --panels urda,cutc          # report a subset of panels\n"
            "  openbiota run --compact                   # one-screen answer\n"
            "  openbiota panels                          # list panels\n"
            "  openbiota doctor                          # check dependencies and hardware\n"
        ),
    )
    parser.add_argument("--version", action="version", version=f"openbiota {__version__}")
    sub = parser.add_subparsers(dest="command")

    run = sub.add_parser("run", help="run the screen (default command)")
    _add_run_arguments(run)

    panels = sub.add_parser("panels", help="list or inspect gene panels")
    panels.add_argument("--panels-dir", type=Path, default=REPO_ROOT / "panels")
    panels.add_argument("--show", metavar="NAME", help="print one panel's full definition")
    panels.add_argument("--json", action="store_true", help="emit JSON")

    build = sub.add_parser("build-db", help="fetch references and build the DIAMOND database")
    build.add_argument("--panels-dir", type=Path, default=REPO_ROOT / "panels")
    build.add_argument("--refs-dir", type=Path, default=REPO_ROOT / "refs")
    build.add_argument("--db-panels", default="all", help="comma-separated panel names, or 'all'")
    build.add_argument("--threads", type=int, default=_default_threads())
    build.add_argument("--diamond", default="diamond")
    build.add_argument("--refresh-refs", action="store_true")
    build.add_argument("--rebuild-db", action="store_true")
    build.add_argument("--quiet", action="store_true")

    doctor = sub.add_parser("doctor", help="check dependencies, hardware and inputs")
    doctor.add_argument("--fastq-dir", type=Path, default=REPO_ROOT / "fastq")
    doctor.add_argument("--panels-dir", type=Path, default=REPO_ROOT / "panels")
    doctor.add_argument("--diamond", default="diamond")

    cohort = sub.add_parser(
        "cohort",
        help="build empirical reference ranges from public metagenomes",
        description=(
            "Downloads read prefixes from a public ENA study, screens each sample, and "
            "reduces the results to percentile ranges. Those ranges are what turn a bare "
            "number into 'below / within / above the usual range'."
        ),
    )
    cohort.add_argument("--study", default=DEFAULT_STUDY, help="ENA study accession")
    cohort.add_argument(
        "--samples", type=int, default=100,
        help=(
            "number of runs to use. The band edges are the 5th and 95th "
            "percentiles, and at n=100 each is set by about five samples, so a "
            "smaller cohort makes the extremes - where the report says something "
            "alarming - rest on two or three people."
        ),
    )
    cohort.add_argument(
        "--runs", type=Path, default=None,
        help=(
            "file of run accessions, one per line, to use instead of the first "
            "--samples; a previous reference_ranges.json's 'runs' list rebuilds "
            "that cohort exactly"
        ),
    )
    cohort.add_argument(
        "--reads", type=int, default=600_000, help="read pairs per cohort sample"
    )
    cohort.add_argument("--workers", type=int, default=10, help="parallel downloads")
    cohort.add_argument(
        "--sampling", choices=SAMPLING_MODES, default="random",
        help=(
            "how each cohort sample's reads are chosen: 'random' streams the whole "
            "file and keeps a uniform random subset (correct); 'prefix' takes the "
            "head of the file (cheap, measurably biased)"
        ),
    )
    cohort.add_argument(
        "--seed", type=int, default=DEFAULT_SEED,
        help="seed for --sampling random, so a cohort rebuilds exactly",
    )
    cohort.add_argument(
        "--transport", choices=TRANSPORTS, default="auto",
        help=(
            "where reads are fetched from: 'auto' prefers the SRA toolkit and falls "
            "back to ENA per run; 'ena' refuses the toolkit, which is what you want "
            "when NCBI is throttling and every fastq-dump costs a timeout first; "
            "'sra' is the mirror image"
        ),
    )
    cohort.add_argument(
        "--min-rpob", type=int, default=300,
        help="drop cohort samples with fewer rpoB fragments than this; too few makes that "
             "sample's values imprecise enough to distort the percentiles (default: 300)",
    )
    cohort.add_argument("--panels-dir", type=Path, default=REPO_ROOT / "panels")
    cohort.add_argument("--refs-dir", type=Path, default=REPO_ROOT / "refs")
    cohort.add_argument("--work-dir", type=Path, default=REPO_ROOT / "refs" / "cohort")
    cohort.add_argument(
        "--out", type=Path, default=REPO_ROOT / "refs" / "reference_ranges.json"
    )
    cohort.add_argument("--threads", type=int, default=_default_threads())
    cohort.add_argument("--diamond", default="diamond")
    cohort.add_argument("--keep-fastq", action="store_true", help="do not delete prefixes")
    cohort.add_argument("--quiet", action="store_true")

    validate = sub.add_parser(
        "validate",
        help="measure sensitivity, specificity and quantitative accuracy",
        description=(
            "Builds synthetic communities from reference genomes with known gene content, "
            "simulates reads, screens them, and compares the result against the known truth. "
            "Produces detection sensitivity, false-positive rate and a calibration curve."
        ),
    )
    validate.add_argument("--panels-dir", type=Path, default=REPO_ROOT / "panels")
    validate.add_argument("--refs-dir", type=Path, default=REPO_ROOT / "refs")
    validate.add_argument(
        "--work-dir", type=Path, default=REPO_ROOT / "refs" / "validation"
    )
    validate.add_argument(
        "--out", type=Path, default=REPO_ROOT / "docs" / "validation_results.json"
    )
    validate.add_argument(
        "--reads", type=int, default=4_000_000, help="read pairs per synthetic community"
    )
    validate.add_argument("--threads", type=int, default=_default_threads())
    validate.add_argument("--diamond", default="diamond")
    validate.add_argument("--seed", type=int, default=20260904)
    validate.add_argument("--error-rate", type=float, default=0.002)
    validate.add_argument("--refresh", action="store_true", help="re-download genomes")
    validate.add_argument("--quiet", action="store_true")

    depth = sub.add_parser(
        "depth-check",
        help="test whether the sequencing depth was sufficient",
        description=(
            "Screens the same sample at increasing depths and reports whether the values "
            "have stopped changing. Answers the question 'would deeper sequencing change "
            "this result?' with a measurement rather than an assertion."
        ),
    )
    depth.add_argument("--fastq-dir", type=Path, default=REPO_ROOT / "fastq")
    depth.add_argument("--r1", type=Path)
    depth.add_argument("--r2", type=Path)
    depth.add_argument("--sample", default=None)
    depth.add_argument("--panels-dir", type=Path, default=REPO_ROOT / "panels")
    depth.add_argument("--refs-dir", type=Path, default=REPO_ROOT / "refs")
    depth.add_argument("--work-dir", type=Path, default=REPO_ROOT / "results")
    depth.add_argument(
        "--out", type=Path, default=None,
        help="default: results/<sample>/depth_check.json",
    )
    depth.add_argument(
        "--fractions", default=",".join(str(f) for f in DEFAULT_FRACTIONS),
        help="comma-separated fractions of full depth to screen",
    )
    depth.add_argument("--threads", type=int, default=_default_threads())
    depth.add_argument("--block-size", type=float, default=8.0)
    depth.add_argument("--diamond", default="diamond")
    depth.add_argument("--quiet", action="store_true")

    taxref = sub.add_parser(
        "taxonomic-cohort",
        help="build the species-level reference cohort from curatedMetagenomicData",
        description=(
            "Downloads uniformly processed taxonomic profiles and curated subject metadata "
            "for public stool metagenomes, and assembles a reference distribution for "
            "profile scoring. One snapshot only, so one profiler version."
        ),
    )
    taxref.add_argument("--snapshot", default=DEFAULT_SNAPSHOT)
    taxref.add_argument("--condition", default="control",
                        help="cMD study_condition to select (default: control)")
    taxref.add_argument("--max-samples", type=int, default=None)
    taxref.add_argument("--min-reads", type=int, default=1_000_000)
    taxref.add_argument("--refs-dir", type=Path, default=REPO_ROOT / "refs")
    taxref.add_argument(
        "--out", type=Path, default=REPO_ROOT / "refs" / "taxonomic_cohort.json"
    )
    taxref.add_argument("--stability", action="store_true",
                        help="also compute the cohort-size stability curve")
    taxref.add_argument("--refresh", action="store_true")
    taxref.add_argument("--quiet", action="store_true")

    runall = sub.add_parser(
        "run-all",
        help="screen every sample in the FASTQ directory and compare them",
        description=(
            "Runs the full pipeline on every mate pair found in the FASTQ directory. Each "
            "sample gets its own directory and its own named report, and a side-by-side "
            "comparison across all of them is printed and written to comparison.json."
        ),
    )
    _add_run_arguments(runall)
    runall.add_argument(
        "--only", help="comma-separated sample names to run (default: all found)"
    )
    runall.add_argument(
        "--skip-existing", action="store_true",
        help="skip samples that already have a report",
    )
    runall.add_argument(
        "--continue-on-error", action="store_true", default=True,
        help="carry on when one sample fails (default)",
    )
    runall.add_argument(
        "--stop-on-error", dest="continue_on_error", action="store_false",
        help="abort the batch on the first failure",
    )
    runall.add_argument(
        "--parallel", type=int, default=1, metavar="N",
        help=(
            "screen N samples at once, each with threads/N threads (default 1). The "
            "aligners are multi-threaded but every run also has single-threaded stretches "
            "(fastp, MetaPhlAn's post-processing, tallying, the PDF); overlapping two or "
            "three samples fills those gaps. Each sample's log goes to results/<sample>/run.log."
        ),
    )

    prune = sub.add_parser(
        "prune",
        help="remove superseded reference databases and stale scratch files",
        description=(
            "Each distinct combination of pipeline version, DIAMOND version, reference set "
            "and residue-check configuration gets its own database directory, keyed by a "
            "fingerprint. Old ones accumulate at about 54 MB each. This removes every "
            "database except the one the current configuration would use, and reports what "
            "it would delete before doing so."
        ),
    )
    prune.add_argument("--panels-dir", type=Path, default=REPO_ROOT / "panels")
    prune.add_argument("--refs-dir", type=Path, default=REPO_ROOT / "refs")
    prune.add_argument("--out", type=Path, default=REPO_ROOT / "results")
    prune.add_argument("--diamond", default="diamond")
    prune.add_argument(
        "--alignments", action="store_true",
        help="also delete per-sample alignment TSVs and rarefaction scratch, which are "
             "large and fully regenerable",
    )
    prune.add_argument(
        "--yes", "-y", action="store_true", help="delete without asking"
    )
    prune.add_argument("--quiet", action="store_true")

    host_cmd = sub.add_parser(
        "build-host-index",
        help="download GRCh38 + PhiX and build the Bowtie2 host-filtering index (one-time)",
    )
    host_cmd.add_argument("--refs-dir", type=Path, default=REPO_ROOT / "refs")
    host_cmd.add_argument("--threads", "-p", type=int, default=_default_threads())
    host_cmd.add_argument("--metaphlan", help="path to a MetaPhlAn 3 executable")
    host_cmd.add_argument("--quiet", action="store_true")

    # Biofilm module (spec v08.0 §15). The verbs follow the spec's semantics:
    # `inventory` enumerates inputs before expensive work, `analyze` produces
    # a full structured result even when modules are unresolved, `explain`
    # reads stored evidence without rerunning anything, `compare` keeps
    # profiles separate rather than pooling them, and `datasets audit`
    # reports assay/label/access compatibility for every registered source.
    bf = sub.add_parser(
        "biofilm",
        help="biofilm-related potential: inventory, analyse, explain, compare",
        description=(
            "Two independent axes - concerning and protective - that are never "
            "combined into one number. Research use; not a diagnostic test."
        ),
    )
    bf_sub = bf.add_subparsers(dest="biofilm_cmd", metavar="{inventory,analyze,"
                                                          "explain,compare,validate,"
                                                          "rank-actions,datasets}")
    for name, helptext in (
        ("inventory", "enumerate available reads, caches and reference compatibility"),
        ("analyze", "produce the full structured biofilm result"),
        ("explain", "show the evidence behind a result without rerunning anything"),
        ("rank-actions", "rank matched evidence for an existing result"),
    ):
        p = bf_sub.add_parser(name, help=helptext)
        p.add_argument("--sample", required=True, help="sample ID or results directory")
        p.add_argument("--results-dir", type=Path, default=REPO_ROOT / "results")
        p.add_argument("--cohort", type=Path,
                       default=REPO_ROOT / "refs" / "taxonomic_cohort.json")
        p.add_argument("--output-format", choices=("json", "text"), default="text")
        if name == "explain":
            p.add_argument("--include-candidate-evidence", action="store_true")
        if name == "rank-actions":
            p.add_argument("--include-preclinical", action="store_true", default=True)
            p.add_argument("--include-opposing", action="store_true", default=True)

    bf_cmp = bf_sub.add_parser("compare", help="compare a recipient with donors")
    bf_cmp.add_argument("--recipient", required=True)
    bf_cmp.add_argument("--donors", nargs="+", required=True)
    bf_cmp.add_argument("--results-dir", type=Path, default=REPO_ROOT / "results")
    bf_cmp.add_argument("--cohort", type=Path,
                        default=REPO_ROOT / "refs" / "taxonomic_cohort.json")
    bf_cmp.add_argument("--output-format", choices=("json", "text"), default="text")

    bf_val = bf_sub.add_parser("validate", help="run implementation checks")
    bf_val.add_argument("--registry", default="biofilm-v08")
    bf_val.add_argument("--suite", default="analytical-and-report")
    bf_val.add_argument("--output-format", choices=("json", "text"), default="text")

    bf_ds = bf_sub.add_parser("datasets", help="audit registered public datasets")
    bf_ds.add_argument("action", choices=("audit",), nargs="?", default="audit")
    bf_ds.add_argument("--registry", default="biofilm-v08")
    bf_ds.add_argument("--output-format", choices=("json", "text"), default="text")

    profiles_cmd = sub.add_parser(
        "profiles", help="list or inspect disease similarity profiles"
    )
    profiles_cmd.add_argument("--profiles-dir", type=Path, default=REPO_ROOT / "profiles")
    profiles_cmd.add_argument("--show", metavar="NAME")
    profiles_cmd.add_argument("--json", action="store_true")

    pval = sub.add_parser(
        "validate-profiles",
        help="score profiles against labelled cohorts (AUC + specificity matrix)",
        description=(
            "Runs the profile scoring engine against curatedMetagenomicData cohorts that "
            "carry disease labels. The colorectal cancer AUC is the build gate: if a "
            "known-good signal cannot be recovered on labelled data, no other profile "
            "number means anything. Also produces the cross-profile specificity matrix."
        ),
    )
    pval.add_argument("--profiles-dir", type=Path, default=REPO_ROOT / "profiles")
    pval.add_argument("--refs-dir", type=Path, default=REPO_ROOT / "refs")
    pval.add_argument(
        "--cohort", type=Path, default=REPO_ROOT / "refs" / "taxonomic_cohort.json"
    )
    pval.add_argument(
        "--conditions", default="CRC,IBD,T2D,adenoma",
        help="comma-separated cMD study_condition labels to score against",
    )
    pval.add_argument("--max-per-arm", type=int, default=400)
    pval.add_argument("--snapshot", default=DEFAULT_SNAPSHOT)
    pval.add_argument(
        "--out", type=Path, default=REPO_ROOT / "docs" / "profile_validation.json"
    )
    pval.add_argument("--quiet", action="store_true")

    from openbiota.fmt.cli import add_parser as _add_fmt_parser
    _add_fmt_parser(sub)

    agem = sub.add_parser(
        "age-model",
        help="build, audit or inspect the estimated stool-microbiome chronological age model",
        description=(
            "Clean-room training of the stool-WGS age ensemble — a linear carriage model "
            "(Huang et al. 2020) and a transformer neural network (Myers et al. 2025) — on "
            "curatedMetagenomicData healthy-control stool profiles (spec v04 §9). "
            "`train` reconstructs the one-baseline-per-participant manifest, runs the "
            "leave-one-study-out / leave-one-country-out audits with every transformation "
            "fitted inside the fold, calibrates split-conformal intervals on untouched "
            "studies, learns the OOD thresholds from held-out folds, and freezes the "
            "adult-only bundle. `audit` prints the frozen model card."
        ),
    )
    agem.add_argument("action", choices=("train", "audit", "fetch"))
    agem.add_argument("--refs-dir", type=Path, default=REPO_ROOT / "refs")
    agem.add_argument("--no-pooled-audit", action="store_true",
                      help="skip the all-age pooled reproduction audit (faster)")
    agem.add_argument("--no-benchmarks", action="store_true",
                      help="skip the random-forest / gradient-boosting / SVR benchmarks (needs scikit-learn)")
    agem.add_argument("--quiet", action="store_true")

    # `openbiota mycobiome`: the fungal module (BUILD_SPEC_v08.2 §13).
    myc = sub.add_parser(
        "mycobiome",
        help="gut mycobiome: references prepare, analyze, score",
        description="Fungal sequence signal, supported fungi, strains and the experimental MHS-E1 score. "
                    "Research use; not a diagnostic test.",
    )
    myc_sub = myc.add_subparsers(dest="mycobiome_cmd", metavar="{references,analyze,score}")
    mref = myc_sub.add_parser("references", help="build the fungal reference lock, sketch DB, species index and strain panels")
    mref.add_argument("verb", choices=["prepare", "status"])
    mref.add_argument("--refs-dir", type=Path, default=REPO_ROOT / "refs")
    mref.add_argument("--threads", type=int, default=16)
    mref.add_argument("--no-panels", action="store_true", help="skip strain panel derivation")
    mana = myc_sub.add_parser("analyze", help="run the mycobiome module for one analysed sample")
    mana.add_argument("--sample", required=True)
    mana.add_argument("--r1", type=Path, required=True)
    mana.add_argument("--r2", type=Path)
    mana.add_argument("--output", type=Path, required=True, help="the sample's results directory (holds results.json)")
    mana.add_argument("--refs-dir", type=Path, default=REPO_ROOT / "refs")
    mana.add_argument("--threads", type=int, default=16)
    msc = myc_sub.add_parser("score", help="recompute MHS-E1 from a stored mycobiome record")
    msc.add_argument("--sample-results", type=Path, required=True)
    msc.add_argument("--model", default="MHS-E1")

    return parser


def _add_run_arguments(p: argparse.ArgumentParser) -> None:
    inputs = p.add_argument_group("inputs")
    inputs.add_argument(
        "--fastq-dir", type=Path, default=REPO_ROOT / "fastq",
        help="directory holding the FASTQ pair (default: ./fastq)",
    )
    inputs.add_argument("--r1", type=Path, help="explicit mate 1 path")
    inputs.add_argument("--r2", type=Path, help="explicit mate 2 path")
    inputs.add_argument("--sample", help="sample name (default: inferred from filenames)")
    inputs.add_argument(
        "--no-prefer-gzip", action="store_true",
        help="use the uncompressed FASTQ when both forms exist (default prefers .gz: "
             "DIAMOND reads gzip natively and it halves disk I/O)",
    )

    panels = p.add_argument_group("panels")
    panels.add_argument("--panels-dir", type=Path, default=REPO_ROOT / "panels")
    panels.add_argument(
        "--panels", default="all",
        help="comma-separated panel names to REPORT, or 'all' (default: all)",
    )
    panels.add_argument(
        "--db-panels", default="all",
        help="comma-separated panel names to include in the search DATABASE, or 'all'. "
             "Keeping this at 'all' is recommended: cross-panel best-hit competition is a "
             "specificity mechanism and a wider database costs almost no runtime.",
    )
    panels.add_argument(
        "--min-fragments", type=int,
        help="override every panel's statistical-stability threshold",
    )
    panels.add_argument(
        "--no-residue-check", action="store_true",
        help="skip active-site residue checks (drops two DIAMOND output columns)",
    )

    perf = p.add_argument_group("performance")
    perf.add_argument("--threads", "-p", type=int, default=_default_threads())
    perf.add_argument(
        "--block-size", "-b", type=float, default=8.0,
        help="DIAMOND query block size in billions of letters; RAM is roughly 6x this in GB "
             "(default: 8)",
    )
    perf.add_argument("--index-chunks", "-c", type=int, default=1)
    perf.add_argument("--sensitivity", choices=SENSITIVITIES, default="default",
                      help="DIAMOND sensitivity; --very-sensitive is ~10x slower and "
                           "unnecessary at 50-60%% identity (default: default)")
    perf.add_argument("--evalue", type=float, default=1e-5)
    perf.add_argument("--diamond", default="diamond", help="path to the diamond executable")
    perf.add_argument(
        "--subsample", type=int, metavar="N",
        help="run against the first N read pairs only (fast smoke test)",
    )

    cache = p.add_argument_group("caching")
    cache.add_argument("--refs-dir", type=Path, default=REPO_ROOT / "refs")
    cache.add_argument("--out", type=Path, default=REPO_ROOT / "results")
    cache.add_argument("--refresh-refs", action="store_true", help="re-download reference sets")
    cache.add_argument("--rebuild-db", action="store_true", help="rebuild the DIAMOND database")
    cache.add_argument("--force", action="store_true", help="ignore cached DIAMOND hit output")

    ext = p.add_argument_group("report extension (v0.8.3)")
    ext.add_argument(
        "--simulate", choices=("auto", "always", "never"), default="auto",
        help="model-assisted synbiotic scenarios: 'auto' solves them (and caches every "
             "linear programme, so a regenerated report is seconds), 'never' reports only "
             "readiness, 'always' forces the solve even when a cache would serve it",
    )
    ext.add_argument(
        "--lp-workers", type=int, default=None,
        help="worker processes for the scenario solver (default: cores minus two, at most 16)",
    )

    qc = p.add_argument_group("input profiling")
    qc.add_argument("--no-qc", action="store_true", help="skip input validation and read stats")
    qc.add_argument(
        "--no-read-count", action="store_true",
        help="skip the exact record count (a full pass over the input)",
    )
    qc.add_argument("--qc-sample-reads", type=int, default=DEFAULT_SAMPLE_READS)
    qc.add_argument("--no-taxonomy", action="store_true",
                    help="skip the rpoB community profile")
    qc.add_argument("--no-fastp", action="store_true",
                    help="skip the fastp pass (duplication, Q30, adapters, insert size)")
    qc.add_argument("--fastp", help="path to a fastp executable")
    qc.add_argument("--trim", action="store_true",
                    help="write adapter/polyG-trimmed reads with fastp and search those "
                         "instead of the originals (changes every downstream cache)")
    qc.add_argument("--detect-adapters", action="store_true",
                    help="run fastp's adapter auto-detection pre-pass (slow; overlap "
                         "trimming already handles paired-end adapters)")

    out = p.add_argument_group("output")
    out.add_argument("--compact", action="store_true", help="print the one-screen answer only")
    out.add_argument("--json", action="store_true", help="print results.json to stdout")
    out.add_argument("--no-panel-hits", action="store_true",
                     help="skip writing per-panel hits_<panel>_<mate>.tsv views")
    out.add_argument("--no-pdf", action="store_true", help="skip the PDF report")
    out.add_argument(
        "--reference-ranges", type=Path, default=REPO_ROOT / "refs" / "reference_ranges.json",
        help="percentile ranges from `openbiota cohort` (default: refs/reference_ranges.json)",
    )
    out.add_argument(
        "--validation", type=Path, default=REPO_ROOT / "docs" / "validation_results.json",
        help="measured performance from `openbiota validate`",
    )
    out.add_argument(
        "--profile-validation", type=Path,
        default=REPO_ROOT / "docs" / "profile_validation.json",
        help="labelled-cohort AUCs and specificity matrix from `openbiota validate-profiles`",
    )
    out.add_argument("--quiet", "-q", action="store_true", help="suppress progress on stderr")

    sim = p.add_argument_group(
        "profile similarity",
        "Taxonomic engine and disease-similarity profiles. Resemblance to a published "
        "group-level pattern; never a diagnosis or a probability.",
    )
    sim.add_argument("--no-profiles", action="store_true",
                     help="skip the taxonomic engine and profile scoring entirely")
    sim.add_argument("--profiles-dir", type=Path, default=REPO_ROOT / "profiles")
    sim.add_argument("--profile", dest="profile_names", default="all",
                     help="comma-separated profile names to score, or 'all'")
    sim.add_argument("--metaphlan", help="path to a MetaPhlAn 3 executable")
    sim.add_argument("--metaphlan4", help="path to a MetaPhlAn 4 executable (extended catalogue)")
    sim.add_argument("--no-assembly", action="store_true",
        help="skip targeted assembly even when SingleM triggers it (for quick reruns; the completeness gate records the skip)")
    sim.add_argument("--lanes", nargs="*", default=None,
        help="run only these expanded lanes (metaphlan_jan26 sylph_globdb motus4 kraken_uhgg singlem_globdb); default all")
    sim.add_argument("--no-extended-catalogue", action="store_true",
                     help="skip the MetaPhlAn 4 SGB inventory (reported, never scored)")
    sim.add_argument("--no-pathogens", action="store_true",
                     help="skip the pathogen screen (bacteria, protozoa, worms, "
                          "microsporidia, fungi, viruses and their determinants)")
    sim.add_argument(
        "--strict-cdiff", action="store_true",
        help="report a Clostridioides difficile species detection as a positive pathogen "
             "finding even when the toxin genes (tcdA/tcdB) are not found, instead of "
             "filing it as carriage. The report states that this mode was used.",
    )
    sim.add_argument("--skip-host-filter", action="store_true",
                     help="profile without removing host reads (host fraction unreported)")
    sim.add_argument(
        "--taxonomic-cohort", type=Path,
        default=REPO_ROOT / "refs" / "taxonomic_cohort.json",
        help="species reference cohort from `openbiota taxonomic-cohort`",
    )
    sim.add_argument(
        "--functional-cohort", type=Path,
        default=REPO_ROOT / "refs" / "reference_cohort_samples.json",
        help="per-sample gene-capacity values from `openbiota cohort`",
    )
    sim.add_argument(
        "--match", default="country,age_band,sex",
        help="reference-matching variables (default: country,age_band,sex); a stratum "
             "below 30 samples falls back to the full cohort and is flagged",
    )
    sim.add_argument("--subject-country", help="ISO3 country code of the subject, e.g. USA")
    sim.add_argument("--subject-age", type=float, help="subject age in years")
    sim.add_argument("--subject-sex", choices=("male", "female"))
    sim.add_argument("--illness-duration-years", type=float,
                     help="for duration-dependent profiles; both strata reported if omitted")
    sim.add_argument("--antibiotics-days-ago", type=int,
                     help="days since the last course of antibiotics, if any")
    sim.add_argument("--onset-date", help="symptom onset date (YYYY-MM-DD), for abstention rules")
    sim.add_argument("--sampling-date", help="stool sampling date (YYYY-MM-DD)")
    sim.add_argument("--stool-form", help="Bristol stool form or a free-text description")
    sim.add_argument("--medications", help="current medications, free text")
    sim.add_argument("--include-opt-in", action="store_true",
                     help="also score opt-in profiles (psychiatric and other sensitive "
                          "conditions that are withheld unless asked for)")
    sim.add_argument("--taxa-dir", type=Path, default=REPO_ROOT / "taxa",
                     help="curated taxon groups (butyrate producers, oral-origin …) for the "
                          "community overview and carrier features")
    sim.add_argument("--mode", choices=("participant", "clinician", "research"), default="research",
                     help="report mode (spec v04 §3). This tool is released to internal research mode "
                          "(the default): the age estimate is produced without an actual age (no residual), "
                          "and evidence cards show studied doses for supplements and diets. participant "
                          "mode fails closed: no age estimate without a known adult age, no doses; "
                          "clinician mode additionally shows prescription/FMT protocols")
    sim.add_argument("--no-age", action="store_true",
                     help="skip the estimated stool-microbiome chronological age (research use)")


def _split_list(value: str | None) -> list[str] | None:
    if not value or value.strip().lower() == "all":
        return None
    return [item.strip() for item in value.split(",") if item.strip()]


# --------------------------------------------------------------------------- #
# commands
# --------------------------------------------------------------------------- #


def cmd_panels(args: argparse.Namespace) -> int:
    panel_set = load_panel_set(args.panels_dir)
    if args.show:
        panel = panel_set.by_name(args.show)
        if args.json:
            print(json.dumps(_panel_json(panel), indent=2))
            return 0
        print(f"panel: {panel.name}{'  [extension]' if panel.extension else ''}")
        print(f"  description: {' '.join(panel.description.split())}")
        print(f"  metabolite:  {panel.metabolite}")
        if panel.pathway:
            print(f"  pathway:     {' '.join(panel.pathway.split())}")
        print(f"  aggregate:   {panel.aggregate} ({' '.join(panel.aggregate_reason.split())})")
        print(f"  stability:   >= {panel.min_fragments_for_stability} fragments")
        print(f"  min aln:     {panel.min_alignment_aa} aa")
        for role, entries in (("targets", panel.targets), ("decoys", panel.decoys)):
            print(f"  {role}:")
            for entry in entries:
                print(f"    - {entry.id} ({entry.label}) min_identity={entry.min_identity}")
                print(f"      query: {entry.query}")
        if panel.residue_check:
            rc = panel.residue_check
            print(
                f"  residue check: {rc.anchor_accession} position {rc.canonical_position} "
                f"in {sorted(rc.accepted_residues)} applied to {list(rc.applies_to)}"
            )
        if panel.upper_bound_reason:
            print(f"  upper bound: {' '.join(panel.upper_bound_reason.split())}")
        for caveat in panel.caveats:
            print(f"  caveat: {' '.join(caveat.split())}")
        print(f"  citation: {' '.join(panel.citation.split())}")
        return 0

    if args.json:
        print(json.dumps({
            "normalizer": {
                "id": panel_set.normalizer.id,
                "query": panel_set.normalizer.query,
                "description": " ".join(panel_set.normalizer_description.split()),
            },
            "panels": [_panel_json(p) for p in panel_set.panels],
        }, indent=2))
        return 0

    print(f"normaliser: {panel_set.normalizer.id} — {' '.join(panel_set.normalizer_description.split())[:70]}")
    print()
    print(f"{'panel':<14}{'metabolite':<44}{'genes':<26}{'agg':<8}flags")
    print("-" * 100)
    for panel in panel_set.panels:
        genes = ",".join(t.gene or t.id for t in panel.targets)
        flags = []
        if panel.extension:
            flags.append("extension")
        if panel.residue_check and panel.residue_check.enabled:
            flags.append("residue-check")
        if panel.upper_bound_reason:
            flags.append("upper-bound")
        print(
            f"{panel.name:<14}{panel.metabolite[:43]:<44}{genes[:25]:<26}"
            f"{panel.aggregate:<8}{','.join(flags)}"
        )
    print()
    print("Deliberately excluded (near-universal across bacteria, so uninformative):")
    print("  ldh (lactate), gadB (GABA), speE (spermidine), murI (D-glutamate).")
    print("Not screenable: urolithin A — the responsible gut bacteria are unidentified.")
    return 0


def _panel_json(panel: Panel) -> dict[str, Any]:
    return {
        "name": panel.name,
        "description": " ".join(panel.description.split()),
        "metabolite": panel.metabolite,
        "pathway": panel.pathway,
        "extension": panel.extension,
        "aggregate": panel.aggregate,
        "min_fragments_for_stability": panel.min_fragments_for_stability,
        "min_alignment_aa": panel.min_alignment_aa,
        "known_producers": list(panel.organisms),
        "targets": [
            {
                "id": e.id, "label": e.label, "gene": e.gene, "query": e.query,
                "min_identity": e.min_identity, "max_sequences": e.max_sequences,
            }
            for e in panel.targets
        ],
        "decoys": [
            {"id": e.id, "label": e.label, "query": e.query, "min_identity": e.min_identity}
            for e in panel.decoys
        ],
        "residue_check": (
            None
            if panel.residue_check is None
            else {
                "enabled": panel.residue_check.enabled,
                "anchor_accession": panel.residue_check.anchor_accession,
                "canonical_position": panel.residue_check.canonical_position,
                "accepted_residues": sorted(panel.residue_check.accepted_residues),
                "applies_to": list(panel.residue_check.applies_to),
            }
        ),
        "upper_bound_reason": panel.upper_bound_reason,
        "caveats": [" ".join(c.split()) for c in panel.caveats],
        "citation": " ".join(panel.citation.split()),
    }


def cmd_doctor(args: argparse.Namespace) -> int:
    ok = True
    print(f"openbiota {__version__}")
    print(f"python              {sys.version.split()[0]} ({sys.executable})")
    # Kept deliberately despite requires-python: `doctor` exists to tell the
    # user what is wrong with their environment, and pip can be persuaded to
    # install into an unsupported interpreter.
    if sys.version_info < (3, 11):  # noqa: UP036
        print("  FAIL: Python 3.11 or newer is required")
        ok = False

    info = platform_info()
    print(f"platform            {info.system}/{info.machine}, {info.cpu_brand}")
    print(f"cpus                {info.logical_cpus} logical")
    print(f"memory              {human_bytes(info.memory_bytes)}")
    if info.translated:
        print("  WARN: this Python is running translated under Rosetta (~2x slowdown)")

    exe = shutil.which(args.diamond)
    if exe is None:
        print(f"diamond             NOT FOUND ({args.diamond!r} not on PATH)")
        print("  install: brew install diamond   (or the brewsci/bio tap, or build from source)")
        print("  avoid bioconda osx-64 on Apple Silicon: it runs under Rosetta and costs ~2x")
        ok = False
    else:
        try:
            version = diamond_version(args.diamond)
        except OpenBiotaError as exc:
            print(f"diamond             BROKEN: {exc}")
            ok = False
        else:
            from openbiota.search import diamond_binary_arch

            arch = diamond_binary_arch(args.diamond)
            print(f"diamond             {version.splitlines()[0]}")
            print(f"                    {exe} [{arch}]")
            if info.machine == "arm64" and arch == "x86_64":
                print("  WARN: x86_64 DIAMOND on arm64 host — runs under Rosetta, ~2x slower")

    try:
        import yaml  # noqa: F401
    except ImportError:
        print("PyYAML              NOT INSTALLED  (pip install PyYAML)")
        ok = False
    else:
        print(f"PyYAML              {yaml.__version__}")

    try:
        panel_set = load_panel_set(args.panels_dir)
    except OpenBiotaError as exc:
        print(f"panels              INVALID: {exc}")
        ok = False
    else:
        n_entries = len(panel_set.all_entries())
        print(f"panels              {len(panel_set.panels)} panels, {n_entries} reference entries")

    try:
        sample_input = discover_sample(fastq_dir=args.fastq_dir)
    except OpenBiotaError as exc:
        print(f"inputs              {exc}")
    else:
        print(f"inputs              sample {sample_input.sample!r}")
        for mate in sample_input.mates:
            print(
                f"                    {mate.label}: {mate.path.name} "
                f"({human_bytes(mate.path.stat().st_size)}{', gzip' if mate.gzipped else ''})"
            )

    print()
    print("OK" if ok else "PROBLEMS FOUND — see above")
    return 0 if ok else 1


def cmd_build_db(args: argparse.Namespace) -> int:
    reporter = Reporter(verbose=not args.quiet)
    panel_set = load_panel_set(args.panels_dir)
    entries = _entries_for_db(panel_set, _split_list(args.db_panels))
    with reporter.stage("reference database"):
        database = build_reference_database(
            panel_set,
            entries,
            refs_dir=args.refs_dir,
            reporter=reporter,
            refresh=args.refresh_refs,
            rebuild=args.rebuild_db,
            diamond=args.diamond,
            threads=args.threads,
        )
    print(f"database:    {database.dmnd_path}")
    print(f"fingerprint: {database.fingerprint}")
    print(f"proteins:    {database.n_sequences:,} ({database.total_residues():,} aa)")
    print()
    print(f"{'entry':<22}{'role':<12}{'seqs':>7}{'mean aa':>9}{'avail':>8}  length window")
    for key, stats in sorted(database.entries.items()):
        print(
            f"{key:<22}{stats.role:<12}{stats.n_sequences:>7,}{stats.mean_length_aa:>9.0f}"
            f"{(stats.total_available if stats.total_available >= 0 else '?'):>8}  "
            f"{stats.length_window}"
        )
    return 0


def _entries_for_db(panel_set: PanelSet, names: Sequence[str] | None) -> list[Any]:
    selected = panel_set.select(names)
    entries = [panel_set.normalizer]
    for panel in selected:
        entries.extend(panel.entries)
    return entries


def screen_one(
    *,
    database: Any,
    panel_set: PanelSet,
    sample_input: Any,
    out_dir: Path,
    config: SearchConfig,
    fields: Sequence[str],
    reporter: Reporter,
    subsample: int | None = None,
    force: bool = False,
) -> Any:
    """Search one sample and tally it. Shared by `run`, `cohort` and `validate`."""
    searches = run_search(
        database=database,
        sample_input=sample_input,
        out_dir=out_dir,
        config=config,
        fields=fields,
        reporter=reporter,
        subsample=subsample,
        force=force,
    )
    return tally(
        streams=[(s.mate, iter_lines(s.hits_path)) for s in searches],
        database=database,
        panel_set=panel_set,
        fields=fields,
    )


def cmd_cohort(args: argparse.Namespace) -> int:
    reporter = Reporter(verbose=not args.quiet)
    panel_set = load_panel_set(args.panels_dir)
    entries = _entries_for_db(panel_set, None)

    with reporter.stage("reference database"):
        database = build_reference_database(
            panel_set,
            entries,
            refs_dir=args.refs_dir,
            reporter=reporter,
            diamond=args.diamond,
            threads=args.threads,
        )

    with reporter.stage("ENA discovery"):
        runs = list_runs(args.study, reporter)
        title = study_title(args.study)
        if title:
            reporter.info(f"  study title: {title}")
    if args.runs is not None:
        wanted = [
            line.strip() for line in Path(args.runs).read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.startswith("#")
        ]
        by_run = {r.run: r for r in runs}
        missing = [w for w in wanted if w not in by_run]
        if missing:
            raise OpenBiotaError(f"{len(missing)} runs in {args.runs} are not in {args.study}: {missing[:5]}")
        chosen = [by_run[w] for w in wanted]
        reporter.info(f"using {len(chosen)} runs listed in {args.runs} at {args.reads:,} read pairs each")
    else:
        chosen = runs[: args.samples]
        reporter.info(f"using the first {len(chosen)} runs by accession at {args.reads:,} read pairs each")

    # Sampling the whole of every run is far too much data to hold at once, so
    # download and screening are overlapped: a bounded pool keeps a few runs
    # ready while the screener works through them, and each run's reads are
    # deleted the moment its numbers are in. Peak disk is the pool, not the
    # cohort.
    fastq_dir = args.work_dir / "fastq"

    fields = output_fields(with_residues=bool(database.anchors))
    config = SearchConfig(
        diamond=args.diamond, threads=args.threads, block_size=4.0, index_chunks=1
    )

    panel_values: dict[str, list[float]] = {p.name: [] for p in panel_set.panels}
    gene_values: dict[str, list[float]] = {
        e.key: [] for p in panel_set.panels for e in p.targets
    }
    per_sample: list[dict[str, Any]] = []

    def acquire(run: CohortRun) -> tuple[CohortRun, Path, Path]:
        return download_cohort(
            [run],
            dest_dir=fastq_dir,
            read_pairs=args.reads,
            workers=1,
            reporter=Reporter(verbose=False),
            sampling=args.sampling,
            seed=args.seed,
            transport=args.transport,
        )[0]

    with (
        reporter.stage(f"sampling and screening {len(chosen)} cohort runs"),
        concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool,
    ):
        reporter.info(
            f"  {args.sampling} sampling at {args.reads:,} read pairs, "
            f"{args.workers} concurrent downloads via {args.transport}"
            + (f", seed {args.seed}" if args.sampling == "random" else "")
        )
        pending = {pool.submit(acquire, run): run for run in chosen}
        index = 0
        for future in concurrent.futures.as_completed(pending):
            index += 1
            queued = pending[future]
            try:
                run, r1, r2 = future.result()
            except (OpenBiotaError, OSError) as exc:
                reporter.warn(f"  [{index}/{len(chosen)}] {queued.run} not sampled — {exc}")
                continue
            sample_input = discover_sample(r1=r1, r2=r2, sample=run.run)
            out_dir = args.work_dir / "results" / run.run
            try:
                result = screen_one(
                    database=database,
                    panel_set=panel_set,
                    sample_input=sample_input,
                    out_dir=out_dir,
                    config=config,
                    fields=fields,
                    reporter=Reporter(verbose=False),
                )
            except OpenBiotaError as exc:
                reporter.warn(f"  {run.run} failed: {exc}")
                continue

            if result.normalizer.fragments < max(1, args.min_rpob):
                reporter.warn(
                    f"  {run.run} excluded: {result.normalizer.fragments} rpoB fragments, "
                    f"below the --min-rpob floor of {args.min_rpob}"
                )
                continue

            row: dict[str, Any] = {
                "run": run.run,
                "sample": run.sample,
                "rpob_fragments": result.normalizer.fragments,
            }
            for panel in result.panels:
                if panel.copies_per_100_genomes is not None:
                    panel_values[panel.panel.name].append(panel.copies_per_100_genomes)
                    row[panel.panel.name] = round(panel.copies_per_100_genomes, 4)
                for gene in panel.targets:
                    if gene.copies_per_100_genomes is not None:
                        gene_values[gene.key].append(gene.copies_per_100_genomes)
            per_sample.append(row)

            reporter.info(
                f"  [{index}/{len(chosen)}] {run.run}: "
                f"rpoB {result.normalizer.fragments:,}, "
                + ", ".join(
                    f"{p.panel.name} {p.copies_per_100_genomes:.1f}"
                    for p in result.panels
                    if p.copies_per_100_genomes is not None
                )[:110]
            )

            if not args.keep_fastq:
                for path in (r1, r2):
                    path.unlink(missing_ok=True)
                    path.with_suffix(path.suffix + ".done").unlink(missing_ok=True)

    usable = len(per_sample)
    if usable < 5:
        raise OpenBiotaError(
            f"only {usable} cohort samples produced usable results; need at least 5 to "
            "build percentile ranges. Re-run with more --samples."
        )

    ranges = build_ranges(
        study=args.study,
        description=(
            f"{usable} adult stool shotgun metagenomes from ENA study {args.study}, "
            + (
                f"a uniform random {args.reads:,} read pairs of each (seed {args.seed})"
                if args.sampling == "random"
                else f"the first {args.reads:,} read pairs of each"
            )
            + ". Same body site, library strategy and NovaSeq X platform family as "
            "the screened sample."
            + (f" ENA describes the study as: {title}." if title else "")
        ),
        study_title=title,
        reads_per_sample=args.reads,
        sampling=args.sampling,
        seed=args.seed if args.sampling == "random" else None,
        runs=[row["run"] for row in per_sample],
        panel_values={k: v for k, v in panel_values.items() if v},
        gene_values={k: v for k, v in gene_values.items() if v},
    )
    ranges.save(args.out)
    (args.out.parent / "reference_cohort_samples.json").write_text(
        json.dumps(per_sample, indent=2) + "\n", encoding="utf-8"
    )

    reporter.ok(f"reference ranges written to {args.out} from {usable} samples")
    print()
    print(f"Reference ranges — {usable} samples from {args.study}")
    print(f"{'panel':<14}{'p5':>9}{'p25':>9}{'median':>9}{'p75':>9}{'p95':>9}{'detected':>10}")
    print("-" * 69)
    for name, pr in sorted(ranges.panels.items()):
        print(
            f"{name:<14}{pr.percentiles[5]:>9.2f}{pr.percentiles[25]:>9.2f}"
            f"{pr.percentiles[50]:>9.2f}{pr.percentiles[75]:>9.2f}{pr.percentiles[95]:>9.2f}"
            f"{pr.detection_rate:>9.0%}"
        )
    return 0


def cmd_run_all(args: argparse.Namespace) -> int:
    """Screen every sample in the FASTQ directory, then compare them."""
    reporter = Reporter(verbose=not args.quiet)
    samples = discover_samples(Path(args.fastq_dir), prefer_gzip=not args.no_prefer_gzip)

    wanted = _split_list(args.only)
    if wanted:
        known = {s.sample for s in samples}
        missing = [n for n in wanted if n not in known]
        if missing:
            raise OpenBiotaError(
                f"sample(s) not found in {args.fastq_dir}: {', '.join(missing)}. "
                f"Available: {', '.join(sorted(known))}"
            )
        samples = [s for s in samples if s.sample in set(wanted)]

    if args.skip_existing:
        remaining = []
        for sample in samples:
            report = Path(args.out) / sample.sample / f"{sample.sample}_report.pdf"
            if report.is_file():
                reporter.record(f"    {sample.sample}: report exists, skipping")
            else:
                remaining.append(sample)
        samples = remaining

    if not samples:
        reporter.warn("nothing to do")
        return 0

    total_bytes = sum(m.path.stat().st_size for s in samples for m in s.mates)
    reporter.info(
        f"batch: {len(samples)} sample(s), {human_bytes(total_bytes)} of input — "
        + ", ".join(s.sample for s in samples)
    )

    outcomes: list[SampleOutcome] = []
    batch_started = time.monotonic()
    if args.parallel > 1 and len(samples) > 1:
        outcomes = _run_all_parallel(args, samples, reporter)
    else:
        for index, sample in enumerate(samples, start=1):
            print()
            print("=" * 88)
            print(f" [{index}/{len(samples)}]  {sample.sample}")
            print("=" * 88)
            started = time.monotonic()
            out_dir = Path(args.out) / sample.sample
            try:
                code = cmd_run(args, sample_input=sample, print_summary=False)
                elapsed = time.monotonic() - started
                results_path = out_dir / "results.json"
                if code != 0 and not results_path.is_file():
                    raise OpenBiotaError(f"run exited with status {code}")
                results = json.loads(results_path.read_text(encoding="utf-8"))
                outcomes.append(
                    outcome_from_results(sample.sample, out_dir, elapsed, results)
                )
                reporter.ok(
                    f"[{index}/{len(samples)}] {sample.sample} in "
                    f"{human_duration(elapsed)} — {out_dir / (sample.sample + '_report.pdf')}"
                )
            except Exception as exc:  # noqa: BLE001 — one bad sample must not lose the batch
                elapsed = time.monotonic() - started
                reporter.warn(f"[{index}/{len(samples)}] {sample.sample} FAILED: {exc}")
                outcomes.append(
                    SampleOutcome(
                        sample=sample.sample, ok=False, out_dir=out_dir,
                        elapsed_s=elapsed, error=f"{type(exc).__name__}: {exc}",
                    )
                )
                if not args.continue_on_error:
                    break

    comparison = write_comparison(Path(args.out) / "comparison.json", outcomes)
    total_elapsed = time.monotonic() - batch_started

    print()
    print("=" * 88)
    print(f" BATCH COMPARISON — {len([o for o in outcomes if o.ok])} of {len(outcomes)} "
          f"samples in {human_duration(total_elapsed)}")
    print("=" * 88)
    print_comparison(outcomes)
    print("  Percentiles are against a reference cohort, not against each other.")
    print("  Disease-pattern percentiles are resemblance to a published group pattern —")
    print("  not a diagnosis and not a probability of disease.")
    print()
    reporter.ok(f"comparison written to {comparison}")
    return 0 if all(o.ok for o in outcomes) else 1


#: run-all options that must not be forwarded to the per-sample ``openbiota run``.
_BATCH_ONLY_FLAGS: Final = {"--only", "--parallel"}
_BATCH_ONLY_SWITCHES: Final = {"--skip-existing", "--continue-on-error", "--stop-on-error"}


def _single_run_argv(sample: SampleInput, threads: int) -> list[str]:
    """Rebuild the command line for one sample from the batch invocation.

    The batch's own options are dropped, the pair is named explicitly so
    discovery is not re-run in a multi-sample directory, and the thread
    budget is the batch's share.
    """
    raw = list(sys.argv[1:])
    out: list[str] = []
    skip_next = False
    for token in raw:
        if skip_next:
            skip_next = False
            continue
        if token == "run-all":
            continue
        if token in _BATCH_ONLY_SWITCHES:
            continue
        if token in _BATCH_ONLY_FLAGS or token in {"--threads", "-p", "--fastq-dir", "--sample", "--r1", "--r2"}:
            skip_next = True
            continue
        if any(token.startswith(f"{flag}=") for flag in (*_BATCH_ONLY_FLAGS, "--threads", "--fastq-dir", "--sample", "--r1", "--r2")):
            continue
        out.append(token)
    argv = ["run", *out, "--sample", sample.sample, "--threads", str(threads), "--r1", str(sample.mates[0].path)]
    if len(sample.mates) > 1:
        argv += ["--r2", str(sample.mates[1].path)]
    return argv


def _run_all_parallel(
    args: argparse.Namespace, samples: list[SampleInput], reporter: Reporter
) -> list[SampleOutcome]:
    """Screen ``samples`` ``args.parallel`` at a time in child processes.

    Each child is a plain ``openbiota run`` with its share of the threads; its
    console output is kept in ``results/<sample>/run.console.log`` so the
    batch console stays readable; the child writes its own ``run.log``. Memory: DIAMOND holds roughly 6 GB per unit of block size and
    the MetaPhlAn 4 index about 10 GB, so two or three at once suit a
    workstation with 64 GB or more.
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed

    workers = max(1, min(args.parallel, len(samples)))
    per_thread = max(1, args.threads // workers)
    reporter.info(f"running {workers} samples at a time, {per_thread} threads each")

    def one(sample: SampleInput) -> SampleOutcome:
        started = time.monotonic()
        out_dir = Path(args.out) / sample.sample
        out_dir.mkdir(parents=True, exist_ok=True)
        # The child writes its own results/<sample>/run.log — header, every
        # timestamped progress line, and the total at the end. The parent must
        # not open that same file: streaming the child's stdio into it used to
        # overwrite the child's log from offset 0 and leave a file with no
        # timings in it. Console output is captured beside it instead, and is
        # what the error tail shows if the child died before writing its log.
        console_path = out_dir / "run.console.log"
        argv = _single_run_argv(sample, per_thread)
        with console_path.open("w", encoding="utf-8") as console:
            proc = subprocess.run(
                [sys.executable, "-m", "openbiota", *argv], stdout=console, stderr=subprocess.STDOUT,
                check=False, env={**os.environ, "PYTHONUNBUFFERED": "1"},
            )
        elapsed = time.monotonic() - started
        results_path = out_dir / "results.json"
        if proc.returncode != 0 and not results_path.is_file():
            tail = "\n".join(
                console_path.read_text(encoding="utf-8", errors="replace").splitlines()[-12:]
            )
            return SampleOutcome(
                sample=sample.sample, ok=False, out_dir=out_dir, elapsed_s=elapsed,
                error=f"exit {proc.returncode}; see {console_path}\n{tail}",
            )
        results = json.loads(results_path.read_text(encoding="utf-8"))
        return outcome_from_results(sample.sample, out_dir, elapsed, results)

    outcomes: list[SampleOutcome] = []
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(one, s): s for s in samples}
        for index, future in enumerate(as_completed(futures), start=1):
            sample = futures[future]
            try:
                outcome = future.result()
            except Exception as exc:  # noqa: BLE001 — one bad sample must not lose the batch
                outcome = SampleOutcome(
                    sample=sample.sample, ok=False, out_dir=Path(args.out) / sample.sample,
                    elapsed_s=0.0, error=f"{type(exc).__name__}: {exc}",
                )
            outcomes.append(outcome)
            # Naming a PDF path is a claim that the PDF exists. A sample whose
            # analysis succeeded but whose report failed to render was being
            # announced with a tick and a path to a stale file from an earlier
            # run, which is how a broken layout survived three batches
            # unnoticed. Check the file, and say so when it is not there.
            pdf_path = outcome.out_dir / f"{sample.sample}_report.pdf"
            if outcome.ok and pdf_path.is_file():
                reporter.ok(
                    f"[{index}/{len(samples)}] {sample.sample} in "
                    f"{human_duration(outcome.elapsed_s)} — {pdf_path}"
                )
            elif outcome.ok:
                reporter.warn(
                    f"[{index}/{len(samples)}] {sample.sample} analysed in "
                    f"{human_duration(outcome.elapsed_s)} but NO PDF was written "
                    f"— see {outcome.out_dir / 'run.log'}"
                )
            else:
                reporter.warn(f"[{index}/{len(samples)}] {sample.sample} FAILED: {outcome.error}")
    order = {s.sample: i for i, s in enumerate(samples)}
    outcomes.sort(key=lambda o: order.get(o.sample, 0))
    return outcomes


def cmd_prune(args: argparse.Namespace) -> int:
    """Remove superseded reference databases and, optionally, alignment files."""
    import shutil

    reporter = Reporter(verbose=not args.quiet)
    panel_set = load_panel_set(args.panels_dir)
    entries = _entries_for_db(panel_set, None)
    version = diamond_version(args.diamond)
    current = database_fingerprint(panel_set, entries, version)

    db_root = Path(args.refs_dir) / "db"
    stale = [
        d for d in sorted(db_root.iterdir())
        if d.is_dir() and d.name != current
    ] if db_root.is_dir() else []

    targets: list[tuple[Path, int]] = [(d, _dir_bytes(d)) for d in stale]

    if args.alignments:
        for pattern in ("*/alignments", "*/depth"):
            for sample_dir in sorted(Path(args.out).glob(pattern)):
                size = _dir_bytes(sample_dir)
                if size:
                    targets.append((sample_dir, size))

    if not targets:
        reporter.ok(f"nothing to prune — the current database is {current}")
        return 0

    total = sum(size for _, size in targets)
    print()
    print(f"  current reference database: {current}  (kept)")
    print()
    print(f"  {'to delete':<58}{'size':>12}")
    print("  " + "-" * 70)
    for path, size in targets:
        try:
            shown = path.relative_to(REPO_ROOT)
        except ValueError:
            shown = path
        print(f"  {str(shown)[:57]:<58}{human_bytes(size):>12}")
    print("  " + "-" * 70)
    print(f"  {'total':<58}{human_bytes(total):>12}")
    print()

    if not args.yes:
        print("  Re-run with --yes to delete. Nothing has been removed.")
        return 0

    for path, _ in targets:
        shutil.rmtree(path, ignore_errors=True)
    reporter.ok(f"freed {human_bytes(total)}")
    return 0


def _dir_bytes(path: Path) -> int:
    return sum(p.stat().st_size for p in path.rglob("*") if p.is_file())


def cmd_build_host_index(args: argparse.Namespace) -> int:
    reporter = Reporter(verbose=not args.quiet)
    tools = locate_tools(explicit=args.metaphlan)
    with reporter.stage("host index"):
        ensure_host_index(tools, Path(args.refs_dir), reporter, args.threads)
    return 0


def cmd_profiles(args: argparse.Namespace) -> int:
    profile_set = load_profile_set(args.profiles_dir)
    if args.show:
        profile = profile_set.by_name(args.show)
        if args.json:
            print(json.dumps(profile.to_json(), indent=2))
            return 0
        print(f"profile: {profile.name}  v{profile.version}  [{profile.status}]")
        print(f"  {profile.label}")
        print()
        for line in textwrap.wrap(" ".join(profile.summary.split()), 96):
            print(f"  {line}")
        print()
        if profile.duration_dependent:
            print("  DURATION DEPENDENT — strata are scored separately, never averaged:")
            for s in profile.strata:
                print(f"    {s.name:<14}{s.label}  ({s.criterion})")
            print()
        for name, module in profile.modules.items():
            share = (
                f"{module.weight_in_combined:.0%} of combined"
                if module.weight_in_combined
                else "excluded from combined"
            )
            print(f"  {name} module — {share}, {len(module.features)} features")
            print(f"    {'feature':<38}{'dir':<11}{'w':>5}{'q':>5}{'s':>5}{'c':>5}{'wqsc':>7}")
            for f in module.features:
                print(
                    f"    {f.name[:37]:<38}{f.direction:<11}"
                    f"{f.w:>5.2f}{f.q:>5.2f}{f.s:>5.2f}{f.c:>5.2f}{f.factor:>7.3f}"
                )
            print()
        for check in profile.cross_engine_checks:
            print(f"  cross-engine: {check.taxonomic_feature} vs {check.functional_feature}")
        for rule in profile.abstain_if:
            print(f"  abstain if: {rule.kind} = {rule.value}")
        print()
        for caveat in profile.caveats:
            for line in textwrap.wrap(" ".join(caveat.split()), 94):
                print(f"  ! {line}")
            print()
        return 0

    if args.json:
        print(json.dumps([p.to_json() for p in profile_set.profiles], indent=2))
        return 0

    print(f"{'profile':<14}{'ver':<8}{'status':<20}{'features':>9}  label")
    print("-" * 100)
    for profile in profile_set.profiles:
        n = len(profile.features())
        print(
            f"{profile.name:<14}{profile.version:<8}{profile.status:<20}{n:>9}  "
            f"{profile.label[:44]}"
        )
    print()
    print("Profile similarity is resemblance to a group-level published pattern.")
    print("It is not a diagnosis and not a probability of disease.")
    return 0


def cmd_age_model(args: argparse.Namespace) -> int:
    """Train, audit or fetch data for the chronological-age model (spec v04 §9)."""
    from openbiota import age as agemod

    reporter = Reporter(verbose=not args.quiet)
    refs = Path(args.refs_dir)
    if args.action == "fetch":
        with reporter.stage("fetching curatedMetagenomicData profiles for the age seed studies"):
            agemod.build_manifest(refs / "cmd", reporter)
            paths = agemod.fetch_profiles(refs / "cmd", reporter)
        reporter.ok(f"{len(paths)} study profiles present under {refs / 'cmd' / 'profiles'}")
        return 0
    if args.action == "train":
        bundle = agemod.train(refs, reporter, pooled_audit=not args.no_pooled_audit, benchmarks=not args.no_benchmarks)
        v = bundle.metadata["validation"]
        print(
            f"\n{agemod.MODEL_ID}: adult leave-one-study-out MAE {v['mae_years']:.2f} y "
            f"(kNN baseline {v['knn_baseline_mae_years']:.2f} y), R² {v['r2']:.2f}, "
            f"calibration slope {v['calibration_slope']:.2f}; eligible ages "
            f"{bundle.eligible_age_range[0]:.0f}–{bundle.eligible_age_range[1]:.0f}"
        )
        return 0
    bundle = agemod.load_bundle(refs)
    if bundle is None:
        raise OpenBiotaError(f"no age model bundle at {agemod.bundle_dir(refs)}; run `openbiota age-model train`")
    card = dict(bundle.metadata)
    card.pop("validation_by_study", None)
    print(json.dumps(card, indent=2, default=str))
    return 0


def cmd_validate_profiles(args: argparse.Namespace) -> int:
    """Score profiles against labelled cohorts. The CRC AUC is the build gate."""
    reporter = Reporter(verbose=not args.quiet)
    profile_set = load_profile_set(args.profiles_dir)
    cache = Path(args.refs_dir) / "cmd"

    if not Path(args.cohort).is_file():
        raise OpenBiotaError(
            f"no taxonomic reference cohort at {args.cohort}; run `openbiota taxonomic-cohort` first"
        )
    with reporter.stage("loading reference cohort"):
        healthy = ReferenceCohort.load(Path(args.cohort))
        reporter.info(
            f"  {healthy.n_samples:,} healthy samples, {healthy.n_taxa:,} species, "
            f"snapshot {healthy.snapshot}"
        )

    conditions = [c.strip() for c in args.conditions.split(",") if c.strip()]
    matrix = SpecificityMatrix()
    results: list[dict[str, Any]] = []

    for condition in conditions:
        with reporter.stage(f"assembling labelled cohort — {condition}"):
            try:
                cases = build_cohort(
                    cache_dir=cache,
                    reporter=Reporter(verbose=False),
                    snapshot=args.snapshot,
                    condition=condition,
                    max_samples=args.max_per_arm,
                    min_reads=1_000_000,
                )
            except OpenBiotaError as exc:
                reporter.warn(f"  {condition}: skipped — {exc}")
                continue
            reporter.info(f"  {condition}: {cases.n_samples} labelled samples")

        # Controls come from the healthy cohort, restricted to the studies that
        # contributed cases, so the comparison stays within-study wherever the
        # source cohort provided both arms. Study batch effects are the single
        # biggest source of inflated AUC in pooled microbiome analyses.
        case_studies = {
            cases.metadata[s].study_name for s in cases.sample_ids if s in cases.metadata
        }
        matched_controls = [
            s for s in healthy.sample_ids
            if s in healthy.metadata and healthy.metadata[s].study_name in case_studies
        ]
        controls = (
            healthy.subset(matched_controls)
            if len(matched_controls) >= MIN_PER_ARM
            else healthy
        )
        note = (
            f"controls drawn from the same {len(case_studies)} source studies as the cases"
            if len(matched_controls) >= MIN_PER_ARM
            else "no within-study controls available; the full healthy cohort was used, so "
                 "study batch effects are not controlled"
        )

        for profile in profile_set.profiles:
            if not profile.features(module="taxonomic"):
                continue
            with reporter.stage(f"scoring {profile.name} against {condition}"):
                scores = score_labelled_cohort(
                    profile=profile,
                    case_cohort=cases,
                    control_cohort=controls,
                    condition=condition,
                    reporter=reporter,
                )
                if scores.n_cases < MIN_PER_ARM or scores.n_controls < MIN_PER_ARM:
                    reporter.warn(
                        f"  {profile.name}/{condition}: too few per arm "
                        f"({scores.n_cases}/{scores.n_controls}); skipped"
                    )
                    continue
                result = evaluate_auc(scores, note=note)
            matrix.add(profile.name, condition, result)
            results.append(result.to_json())
            ci = (
                f" [{result.ci_low:.2f}-{result.ci_high:.2f}]"
                if result.ci_low is not None
                else ""
            )
            reporter.info(
                f"  {profile.name} vs {condition}: AUC "
                f"{('n/a' if result.auc is None else f'{result.auc:.3f}')}{ci}, "
                f"study-holdout mean "
                f"{('n/a' if result.holdout_auc is None else f'{result.holdout_auc:.3f}')}"
            )

    payload = {
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "reference_manifest": healthy.manifest(),
        "conditions_scored": conditions,
        "results": results,
        "specificity_matrix": matrix.to_json(),
        "specificity_verdicts": {
            p.name: spread_verdict(matrix, p.name) for p in profile_set.profiles
        },
        "method": (
            "Cases and controls scored through the identical engine. The reference "
            "distribution is built from the healthy cohort only, never from the labelled "
            "data being scored. Leave-one-study-out AUC is reported alongside the pooled "
            "figure because pooling cohorts and splitting randomly leaks recruitment-batch "
            "signal."
        ),
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    print()
    print("=" * 88)
    print(" PROFILE VALIDATION — scoring engine against labelled cohorts")
    print("=" * 88)
    print()
    for entry in results:
        ci = (
            f"  95% CI {entry['auc_ci95'][0]:.2f}-{entry['auc_ci95'][1]:.2f}"
            if entry.get("auc_ci95")
            else ""
        )
        auc = "n/a" if entry["auc"] is None else f"{entry['auc']:.3f}"
        print(
            f"  {entry['profile']:<12} vs {entry['condition']:<10} AUC {auc}"
            f"{ci}   n={entry['n_cases']}/{entry['n_controls']}"
        )
    print()
    print("  CROSS-PROFILE SPECIFICITY MATRIX  (AUC, rows = profile, cols = condition)")
    matrix.print_table()
    print()
    for name, verdict in payload["specificity_verdicts"].items():
        print(f"  {name}: {verdict}")
    print()
    reporter.ok(f"profile validation written to {args.out}")
    return 0


def cmd_taxonomic_cohort(args: argparse.Namespace) -> int:
    """Assemble the species-level reference cohort."""
    reporter = Reporter(verbose=not args.quiet)
    cache = Path(args.refs_dir) / "cmd"

    with reporter.stage("reference cohort (curatedMetagenomicData)"):
        cohort = build_cohort(
            cache_dir=cache,
            reporter=reporter,
            snapshot=args.snapshot,
            condition=args.condition,
            max_samples=args.max_samples,
            min_reads=args.min_reads,
            refresh=args.refresh,
        )
    cohort.save(Path(args.out))

    manifest_path = Path(args.out).with_name(
        Path(args.out).stem + "_manifest.json"
    )
    manifest = cohort.manifest()

    if args.stability:
        with reporter.stage("stability curve"):
            references = build_taxon_references(cohort)
            common = sorted(
                references,
                key=lambda t: -references[t].prevalence,
            )[:40]
            manifest["stability_curve"] = stability_curve(references, taxa=common)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    composition = manifest["composition"]
    print()
    print(f"Taxonomic reference cohort — snapshot {cohort.snapshot} ({cohort.profiler})")
    print(f"  manifest        {cohort.manifest_id}")
    print(f"  samples         {composition['n_samples']:,} from {composition['n_studies']} studies")
    print(f"  species         {cohort.n_taxa:,}")
    print(f"  median depth    {composition['median_reads']:,.0f} reads"
          if composition["median_reads"] else "  median depth    unknown")
    print()
    print(f"  {'country':<24}{'n':>7}")
    for country, count in list(composition["countries"].items())[:10]:
        print(f"  {country:<24}{count:>7,}")
    print()
    print(f"  {'age band':<24}{'n':>7}")
    for band, count in sorted(composition["age_bands"].items()):
        print(f"  {band:<24}{count:>7,}")
    if manifest.get("stability_curve"):
        print()
        print("  stability of the median estimate as the cohort grows")
        print(f"  {'cohort size':<16}{'sd of estimate':>18}")
        for row in manifest["stability_curve"]:
            print(f"  {row['cohort_size']:<16,}{row['median_estimate_sd']:>18.4f}")
    print()
    reporter.ok(f"cohort written to {args.out}, manifest to {manifest_path.name}")
    return 0


def cmd_depth_check(args: argparse.Namespace) -> int:
    """Screen the sample at increasing depths and report convergence."""
    started = time.monotonic()
    reporter = Reporter(verbose=not args.quiet)
    panel_set = load_panel_set(args.panels_dir)

    sample_input = discover_sample(
        fastq_dir=args.fastq_dir if args.r1 is None else None,
        r1=args.r1,
        r2=args.r2,
        sample=args.sample,
        reporter=reporter,
    )

    with reporter.stage("reference database"):
        database = build_reference_database(
            panel_set,
            _entries_for_db(panel_set, None),
            refs_dir=args.refs_dir,
            reporter=reporter,
            diamond=args.diamond,
            threads=args.threads,
        )

    with reporter.stage("counting input"):
        qc = profile_sample_cached(
            [(m.label, m.path) for m in sample_input.mates],
            cache_dir=Path(args.refs_dir) / "cache",
            reporter=reporter,
        )
    full_pairs = qc.read_pairs or 0
    if full_pairs <= 0:
        raise OpenBiotaError("could not determine the input read count")

    fractions = sorted({
        float(f) for f in args.fractions.split(",") if f.strip()
    } | {1.0})
    depths = sorted({max(1, int(full_pairs * f)) for f in fractions})
    reporter.info(
        f"full depth {full_pairs:,} read pairs; screening at "
        + ", ".join(f"{d:,}" for d in depths)
    )

    fields = output_fields(with_residues=bool(database.anchors))
    config = SearchConfig(
        diamond=args.diamond, threads=args.threads, block_size=args.block_size
    )
    points: list[DepthPoint] = []

    for index, depth in enumerate(depths, start=1):
        at_full = depth >= full_pairs
        label = "full depth" if at_full else f"{depth:,} read pairs"
        with reporter.stage(f"depth {index}/{len(depths)} — {label}"):
            # Rarefaction scratch lives under the sample rather than beside it:
            # five sibling `SAMPLE.subsampleN` directories at the top of
            # results/ make the output look like five separate samples.
            out_dir = (
                Path(args.work_dir) / sample_input.sample
                if at_full
                else Path(args.work_dir) / sample_input.sample / "depth" / f"{depth}"
            )
            result = screen_one(
                database=database,
                panel_set=panel_set,
                sample_input=sample_input,
                out_dir=out_dir,
                config=config,
                fields=fields,
                reporter=reporter,
                subsample=None if at_full else depth,
            )
        for panel in result.panels:
            points.append(
                DepthPoint(
                    panel=panel.panel.name,
                    read_pairs=depth,
                    fragments=panel.accepted_fragments,
                    rpob_fragments=result.normalizer.fragments,
                    copies_per_100=panel.copies_per_100_genomes,
                )
            )
        reporter.info(
            f"  rpoB {result.normalizer.fragments:,}; "
            + ", ".join(
                f"{p.panel.name} {p.copies_per_100_genomes:.2f}"
                for p in result.panels
                if p.copies_per_100_genomes is not None
            )[:110]
        )

    report = DepthReport(
        sample=sample_input.sample,
        full_read_pairs=full_pairs,
        fractions=fractions,
        points=points,
        convergence=[
            assess_convergence(points, panel.name) for panel in panel_set.panels
        ],
        generated=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        elapsed_s=time.monotonic() - started,
    )
    out_path = Path(args.out) if args.out else (
        Path(args.work_dir) / sample_input.sample / "depth_check.json"
    )
    report.save(out_path)
    print_depth_summary(report, reporter)
    reporter.ok(f"depth check written to {out_path}")
    return 0


def cmd_validate(args: argparse.Namespace) -> int:
    from openbiota.validate import run_validation

    return run_validation(args)


def _apply_overrides(panel_set: PanelSet, args: argparse.Namespace) -> PanelSet:
    panels = list(panel_set.panels)
    if args.min_fragments is not None:
        panels = [
            dataclasses.replace(p, min_fragments_for_stability=args.min_fragments) for p in panels
        ]
    if args.no_residue_check:
        panels = [
            dataclasses.replace(
                p,
                residue_check=(
                    None
                    if p.residue_check is None
                    else dataclasses.replace(p.residue_check, enabled=False)
                ),
            )
            for p in panels
        ]
    return dataclasses.replace(panel_set, panels=tuple(panels))


def _strain_block(
    out_dir: Any, sample: str, extended: Any, strain_mod: Any
) -> dict[str, Any]:
    """Strain-resolution summary for results.json, or an honest absence.

    The marker lane is a second pass that may not have been run for a given
    sample. When its artifacts are missing this reports `not_run` with the
    command that would produce them — it never reports that the sample has
    no resolvable strains, which is a different statement entirely (§3.2).
    """
    from pathlib import Path as _Path

    hits = sorted(_Path(out_dir).glob("strain/consensus_markers/*.json.bz2"))
    if not hits:
        return {
            "status": "not_run",
            "what_this_is": (
                "Dominant population fingerprints per organism, from marker-gene "
                "consensus sequences. This second pass has not been run for this sample, "
                "so strain resolution is unmeasured rather than absent."
            ),
            "how_to_run": (
                "Align the reads to the marker database emitting SAM, then run "
                "sample2markers.py against the matching 4.1.1 database pkl."
            ),
            "organisms_resolved": 0,
            "organisms_detected_unresolved": 0,
        }

    fingerprints, database = strain_mod.load_fingerprints(hits[0], sample_id=sample)
    species_of = {
        str(row.get("sgb")): str(row.get("species") or "")
        for row in ((extended or {}).get("sgbs") or [])
        if isinstance(row, dict)
    }
    resolved = {k: v for k, v in fingerprints.items() if v.resolved}
    unresolved = {k: v for k, v in fingerprints.items() if not v.resolved}
    rows = [
        {
            "sgb": sgb,
            "species": species_of.get(sgb, ""),
            "markers_resolved": fp.n_markers,
            "callable_bases": fp.callable_bases,
            "median_breadth_percent": round(fp.median_breadth, 1),
            "median_depth": round(fp.median_depth, 2),
            "is_named_strain": False,
        }
        for sgb, fp in sorted(resolved.items(), key=lambda kv: -kv[1].callable_bases)
    ]
    return {
        "status": "resolved",
        "database_release": database,
        "method": "metaphlan 4.1.1 marker consensus (sample2markers)",
        "organisms_resolved": len(resolved),
        "organisms_detected_unresolved": len(unresolved),
        "total_callable_bases": sum(v.callable_bases for v in resolved.values()),
        "organisms": rows,
        "limits": [
            strain_mod.DOMINANT_CONSENSUS_LIMIT,
            strain_mod.MARKER_SCOPE_LIMIT,
        ],
        "what_this_is": (
            "For each organism with enough marker coverage, the dominant population's "
            "consensus across marker genes. This is finer than species and finer than an "
            "SGB label: it is this sample's own population of that organism, which is what "
            "makes donor-recipient comparison possible. It is not a named strain."
        ),
    }


def _resolution_census(sample: str, results_json: dict[str, Any]) -> dict[str, Any]:
    """Every registered target's terminal state, for one sample.

    This is the object that makes an incomplete run impossible to hide. It
    gathers the RA targets, the measured mechanism panels, the non-bacterial
    targets and the uninstalled tool adapters, and reports what each reached.
    """
    from openbiota.resolution import adapters, kingdoms, mechanisms
    from openbiota.resolution import registry as target_registry
    from openbiota.resolution import schema as resolution_schema
    from openbiota.resolution import typing as typing_mod

    calls = list(mechanisms.calls_for_sample(sample, results_json))
    calls += [
        resolution_schema.ResolutionCall(
            sample_id=sample,
            target_id=str(row["target_id"]),
            identity_kind=str(row["identity_kind"]),
            assay_status=str(row["assay_status"]),
            call_id=row.get("call_id"),
            analytical_call=str(row["analytical_call"]),
            reason_codes=tuple(row.get("reason_codes") or ()),
            coverage=resolution_schema.Coverage(**{
                k: v for k, v in (row.get("coverage") or {}).items()
                if k in resolution_schema.Coverage.__slots__
            }),
            plain=str(row.get("plain") or ""),
        )
        for row in ((results_json.get("locus_resolution") or {}).get("targets") or [])
    ]
    calls += kingdoms.all_calls(sample)
    # Every adapter, installed or not: an installed tool nobody ran is a
    # pending measurement and has to be visible as one.
    calls += adapters.calls_for_census(sample)
    # Every RA target, including the operational ones. Skipping those made
    # the census silently smaller than the registry it claims to enumerate:
    # an operational target with no detector run is outstanding work, not
    # a completed screen.
    for target in target_registry.RA_TARGETS:
        operational = target.operational
        calls.append(
            resolution_schema.ResolutionCall(
                sample_id=sample,
                target_id=target.target_id,
                identity_kind=target.identity_kind,
                assay_status=resolution_schema.SCHEDULED,
                reason_codes=(
                    ("target_operational_not_run",) if operational
                    else ("reference_readiness:" + target.reference_readiness,)
                ),
                plain=(
                    (
                        f"{target.label}: the detector is available but was not "
                        "run on this sample, so the target was not measured."
                    )
                    if operational
                    else (target.notes or target.label)
                ),
            )
        )
    for scheme_id in typing_mod.BY_SCHEME:
        calls.append(typing_mod.call_type(
            sample_id=sample, scheme_id=scheme_id, observations=[]
        ))

    census = resolution_schema.census(calls)
    census["targets"] = [c.to_json() for c in calls]
    census["typing_schemes"] = typing_mod.registry_json()
    census["non_bacterial_targets"] = kingdoms.registry_json()["by_capability"]
    census["tool_adapters"] = adapters.registry_json()
    census["ra_targets"] = target_registry.registry_json()
    return census


def cmd_run(
    args: argparse.Namespace,
    *,
    sample_input: SampleInput | None = None,
    print_summary: bool = True,
) -> int:
    """Screen one sample.

    ``sample_input`` lets batch mode supply an already-resolved pair rather
    than re-running discovery, which would fail in a multi-sample directory.
    """
    started = time.monotonic()
    # Wall-clock start for the log header. `time.strftime()` at write time
    # was stamping a ten-minute run as having "started" seven seconds before
    # it finished.
    started_at = time.strftime("%Y-%m-%d %H:%M:%S")
    reporter = Reporter(verbose=not args.quiet)
    timings: dict[str, float] = {}

    reporter.info(f"OpenBiota Gut Health Test (openbiota) v{__version__}")
    check_platform(reporter, args.diamond)

    # ---- panels ---------------------------------------------------------- #
    panel_set = _apply_overrides(load_panel_set(args.panels_dir), args)
    db_entries = _entries_for_db(panel_set, _split_list(args.db_panels))
    # Only panels whose entries are in the search database can be tallied;
    # `--db-panels` narrows the whole run to those (default: every panel).
    panel_set = panel_set.restrict(_split_list(args.db_panels))
    selected_panels = panel_set.select(_split_list(args.panels))
    reporter.info(
        f"panels: {len(selected_panels)} selected for reporting, "
        f"{len(db_entries)} reference entries in the search database"
    )

    # ---- inputs ---------------------------------------------------------- #
    if sample_input is None:
        sample_input = discover_sample(
            fastq_dir=args.fastq_dir if args.r1 is None else None,
            r1=args.r1,
            r2=args.r2,
            sample=args.sample,
            prefer_gzip=not args.no_prefer_gzip,
            reporter=reporter,
        )
    total_input = sum(m.path.stat().st_size for m in sample_input.mates)
    reporter.info(
        f"sample {sample_input.sample!r}: {len(sample_input.mates)} mate file(s), "
        f"{human_bytes(total_input)} on disk"
        + (f" — SUBSAMPLED to the first {args.subsample:,} read pairs" if args.subsample else "")
    )

    out_dir = Path(args.out) / (
        sample_input.sample
        if not args.subsample
        else f"{sample_input.sample}.subsample{args.subsample}"
    )
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- input profiling -------------------------------------------------- #
    qc: SampleQC | None = None
    if not args.no_qc:
        t0 = time.monotonic()
        with reporter.stage("input validation"):
            qc = profile_sample_cached(
                [(m.label, m.path) for m in sample_input.mates],
                cache_dir=Path(args.refs_dir) / "cache",
                sample_reads=args.qc_sample_reads,
                count_all=not args.no_read_count and args.subsample is None,
                reporter=reporter,
                refresh=args.force,
            )
        timings["input validation"] = time.monotonic() - t0
        if args.subsample is not None:
            # The exact record count is skipped for a subsampled run, so fill
            # in the counts that were actually searched rather than leave the
            # totals blank.
            qc = dataclasses.replace(
                qc,
                mates=tuple(
                    dataclasses.replace(
                        m,
                        records=args.subsample,
                        estimated_bases=int(args.subsample * m.mean_length),
                    )
                    for m in qc.mates
                ),
                record_counts_match=True,
                notes=(*qc.notes, f"Subsampled run: first {args.subsample:,} read pairs only."),
            )

    # ---- read preprocessing (fastp) --------------------------------------- #
    fastp_result: FastpResult | None = None
    provenance: HeaderProvenance | None = None
    try:
        provenance = read_header_provenance(sample_input.mates[0].path)
        if provenance.batch_key:
            reporter.record(f"    sequencing batch {provenance.batch_key} ({provenance.instrument})")
    except OSError as exc:  # unreadable header is a QC finding, not a crash
        reporter.warn(f"could not read FASTQ headers: {exc}")
    if not args.no_fastp and args.subsample is None:
        t0 = time.monotonic()
        try:
            with reporter.stage("read preprocessing (fastp)"):
                fastp_tools = locate_fastp(explicit=args.fastp)
                fastp_result = run_fastp(
                    fastp_tools,
                    r1=sample_input.mates[0].path,
                    r2=sample_input.mates[1].path if len(sample_input.mates) > 1 else None,
                    work_dir=out_dir / "preprocess",
                    reporter=reporter,
                    threads=args.threads,
                    trim=args.trim,
                    adapter_detect=args.detect_adapters,
                    force=args.force,
                )
            if args.trim and fastp_result.trimmed:
                sample_input = SampleInput(
                    sample=sample_input.sample,
                    mates=tuple(
                        MateInput(label=m.label, path=p)
                        for m, p in zip(sample_input.mates, fastp_result.trimmed, strict=True)
                    ),
                )
                reporter.info(
                    f"searching trimmed reads: {fastp_result.pairs_out:,} of "
                    f"{fastp_result.pairs_in:,} pairs kept"
                )
        except DependencyError as exc:
            reporter.warn(f"fastp unavailable: {exc}")
        except OpenBiotaError as exc:
            reporter.warn(f"fastp failed: {exc}")
        timings["read preprocessing"] = time.monotonic() - t0

    # ---- reference database ---------------------------------------------- #
    t0 = time.monotonic()
    with reporter.stage("reference database"):
        database = build_reference_database(
            panel_set,
            db_entries,
            refs_dir=args.refs_dir,
            reporter=reporter,
            refresh=args.refresh_refs,
            rebuild=args.rebuild_db,
            diamond=args.diamond,
            threads=args.threads,
        )
    timings["reference database"] = time.monotonic() - t0

    # ---- search ----------------------------------------------------------- #
    needs_residues = any(
        p.residue_check is not None and p.residue_check.enabled for p in panel_set.panels
    ) and bool(database.anchors)
    fields = output_fields(with_residues=needs_residues)
    reporter.info(f"DIAMOND output fields: {' '.join(fields)}")

    config = SearchConfig(
        diamond=args.diamond,
        threads=args.threads,
        block_size=args.block_size,
        index_chunks=args.index_chunks,
        sensitivity=args.sensitivity,
        evalue=args.evalue,
    )
    t0 = time.monotonic()
    with reporter.stage("translated search (diamond blastx)"):
        searches = run_search(
            database=database,
            sample_input=sample_input,
            out_dir=out_dir,
            config=config,
            fields=fields,
            reporter=reporter,
            subsample=args.subsample,
            force=args.force,
        )
    timings["translated search"] = time.monotonic() - t0

    # ---- tally ------------------------------------------------------------ #
    t0 = time.monotonic()
    with reporter.stage("tally and normalisation"):
        result = tally(
            streams=[(s.mate, iter_lines(s.hits_path)) for s in searches],
            database=database,
            panel_set=panel_set,
            fields=fields,
        )
    timings["tally"] = time.monotonic() - t0
    reporter.info(
        f"{result.fragments_with_hit:,} fragments with a hit from "
        f"{result.total_hit_lines:,} read hits "
        f"({result.total_hit_lines - result.fragments_with_hit:,} collapsed as mates); "
        f"rpoB denominator {result.normalizer.fragments:,}"
    )

    # ---- per-panel hit views ---------------------------------------------- #
    if not args.no_panel_hits:
        t0 = time.monotonic()
        with reporter.stage("per-panel hit views"):
            written = write_panel_hit_views(
                hit_files=[(s.mate, s.hits_path) for s in searches],
                database=database,
                panel_names=[p.name for p in selected_panels],
                fields=fields,
                out_dir=out_dir,
            )
        timings["panel hit views"] = time.monotonic() - t0
        reporter.record(f"    wrote {len(written)} per-panel hit view(s)")

    # ---- community profile ------------------------------------------------ #
    profile = None
    if not args.no_taxonomy:
        profile = build_profile(
            phylum_counts=result.normalizer.phylum_counts,
            genus_counts=result.normalizer.genus_counts,
            organism_counts=result.normalizer.organism_counts,
        )

    # ---- taxonomic engine + profile similarity ---------------------------- #
    similarity: SimilarityStage | None = None
    extended: ExtendedCatalogue | None = None
    detection_block: dict[str, Any] = {"schema_version": "openbiota.detection/1.0", "lanes": {}, "lane_status": {},
                                       "reference_releases": []}
    genome_profile: Any = None
    pathogens: PathogenBranchResult | None = None
    age_result: AgeResult | None = None
    age_meta: dict[str, Any] | None = None
    known_total_pairs: int | None = None
    # What is known about the person, what is assumed in its absence, and
    # what stays unknown. One ledger feeds the profile abstention rules, the
    # intervention safety rules and the reading-guide page alike.
    context_ledger = resolve_context(
        mode=args.mode,
        medications=_split_list(args.medications) if args.medications else None,
        antibiotics_days_ago=args.antibiotics_days_ago,
        subject_age=args.subject_age, subject_sex=args.subject_sex, subject_country=args.subject_country,
    )
    if context_ledger.assumptions:
        reporter.record(
            f"    subject context: {len(context_ledger.assumptions)} assumptions in {args.mode} mode "
            f"({', '.join(a.sentence() for a in context_ledger.assumptions[:3])}"
            + (f", +{len(context_ledger.assumptions) - 3} more" if len(context_ledger.assumptions) > 3 else "")
            + "); nothing supplied is withheld, every assumption is printed"
        )

    # The read paths, bound unconditionally. They come straight from the
    # validated input and do not depend on any optional stage, but they used
    # to be assigned inside the MetaPhlAn block: with --no-profiles, or when
    # locate_tools raised, the locus-resolution stage below then hit an
    # UnboundLocalError and the run ended with no results.json, no summary
    # and no run.log, because the reporter buffers until the end.
    r1 = sample_input.mates[0].path
    r2 = sample_input.mates[1].path if len(sample_input.mates) > 1 else None

    if not args.no_profiles:
        t0 = time.monotonic()
        subject = Subject(
            country=args.subject_country,
            age=args.subject_age,
            sex=args.subject_sex,
            illness_duration_years=args.illness_duration_years,
            antibiotics_days_ago=args.antibiotics_days_ago,
            metadata={
                **context_ledger.metadata(),
                "antibiotics": args.antibiotics_days_ago is not None,
                "onset_date": args.onset_date,
                "sampling_date": args.sampling_date,
                "stool_form": args.stool_form,
                "medications": args.medications,
                # The one confounder the sequencing itself can supply: the
                # flowcell:lane batch key read from the FASTQ headers.
                "sequencing_run": (
                    provenance.batch_key
                    if provenance is not None and provenance.batch_key
                    else None
                ),
            },
            assumed_keys=context_ledger.assumed_keys,
        )
        taxonomy = None
        try:
            with reporter.stage("taxonomic engine (MetaPhlAn)"):
                tools = locate_tools(explicit=args.metaphlan)
                # r1/r2 are bound outside this block; see the note there.
                known_total_pairs = (
                    fastp_result.pairs_out if fastp_result is not None and args.trim and fastp_result.trimmed
                    else fastp_result.pairs_in if fastp_result is not None
                    else qc.read_pairs if qc is not None else None
                )
                taxonomy = run_taxonomic_engine(
                    tools=tools, r1=r1, r2=r2, refs_dir=Path(args.refs_dir), out_dir=out_dir,
                    reporter=reporter, threads=args.threads,
                    skip_host_filter=args.skip_host_filter,
                    known_total_pairs=known_total_pairs,
                )
        except DependencyError as exc:
            reporter.warn(f"taxonomic engine unavailable: {exc}")
        except OpenBiotaError as exc:
            reporter.warn(f"taxonomic engine failed: {exc}")

        # Extended catalogue: MetaPhlAn 4 SGBs on the same host-filtered reads.
        # Reported, never scored (see openbiota.engines.metaphlan4).
        if not args.no_extended_catalogue:
            t1 = time.monotonic()
            try:
                with reporter.stage("extended catalogue (MetaPhlAn 4 SGB)"):
                    tools4 = locate_tools4(explicit=args.metaphlan4)
                    if taxonomy is not None and taxonomy.host is not None:
                        q1, q2 = taxonomy.host.nonhost_r1, taxonomy.host.nonhost_r2
                    else:
                        q1 = sample_input.mates[0].path
                        q2 = sample_input.mates[1].path if len(sample_input.mates) > 1 else None
                    # The same alignment is written out as SAM for strain
                    # typing, so the strain lane costs one extra step rather
                    # than a second alignment.
                    strain_dir = out_dir / "strain"
                    marker_sam = strain_dir / f"{sample_input.sample}.markers.sam.bz2"
                    extended = run_metaphlan4(
                        tools4, r1=q1, r2=q2, db_dir=Path(args.refs_dir) / "metaphlan4_db",
                        work_dir=out_dir / "taxonomy", reporter=reporter, threads=args.threads,
                        sam_out=marker_sam,
                    )
            except DependencyError as exc:
                reporter.warn(f"extended catalogue skipped: {exc}")
            except OpenBiotaError as exc:
                reporter.warn(f"extended catalogue failed: {exc}")
            timings["extended catalogue"] = time.monotonic() - t1

            # Strain typing: the marker consensus. This was a second pass run
            # by hand for the first five samples, so the sixth report said
            # "strain typing did not run". It is a stage of every run now,
            # and the completeness gate at the end checks that it produced.
            if extended is not None:
                t1 = time.monotonic()
                try:
                    with reporter.stage("strain typing (marker consensus)"):
                        run_sample2markers(
                            tools4, sam=marker_sam, strain_dir=strain_dir,
                            db_dir=Path(args.refs_dir) / "metaphlan4_db",
                            reporter=reporter, threads=args.threads,
                        )
                except DependencyError as exc:
                    reporter.warn(f"strain typing skipped: {exc}")
                except OpenBiotaError as exc:
                    reporter.warn(f"strain typing failed: {exc}")
                timings["strain typing"] = time.monotonic() - t1

            # Expanded detection (spec 0.8.4): five more lanes over the same
            # host-filtered reads. Each is its own stage and its own record;
            # the inventory merges them, the completeness gate checks them.
            t1 = time.monotonic()
            try:
                from openbiota.expansion import lanes as _lanes

                if taxonomy is not None and taxonomy.host is not None:
                    e1, e2 = taxonomy.host.nonhost_r1, taxonomy.host.nonhost_r2
                else:
                    e1 = sample_input.mates[0].path
                    e2 = sample_input.mates[1].path if len(sample_input.mates) > 1 else None
                detection_block = _lanes.run_all(
                    sample=str(sample_input.sample), r1=e1, r2=e2, refs_dir=Path(args.refs_dir),
                    out_dir=out_dir, threads=args.threads, reporter=reporter,
                    only=tuple(args.lanes) if getattr(args, "lanes", None) else None,
                )
            except Exception as exc:  # noqa: BLE001 - the block records its own failures
                reporter.warn(f"expanded detection unavailable: {type(exc).__name__}: {exc}")
                detection_block = {"schema_version": "openbiota.detection/1.0", "lanes": {}, "lane_status": {},
                                   "reference_releases": [], "error": f"{type(exc).__name__}: {exc}"}
            # The report states measured capability (spec §7): the newest
            # benchmark summary on disk, copied into the results so the PDF
            # reads it from the same file as everything else.
            try:
                from openbiota.expansion import lanes as _lanes_mod

                detection_block["measured_capability"] = _lanes_mod.measured_capability()
            except Exception:  # noqa: BLE001 - absence is recorded as absence
                detection_block["measured_capability"] = {}
            timings["expanded detection"] = time.monotonic() - t1

        # Whole-genome lane: containment profiling of the same reads against
        # every bacterial and archaeal species cluster in GTDB. A detection
        # lane - it resolves organisms the marker catalogues have no markers
        # for - whose abundances stay a secondary reading. Absence of the
        # tool or database records the lane as not assessed, never as zero.
        genome_profile = None
        if not args.no_extended_catalogue:
            t1 = time.monotonic()
            try:
                from openbiota.engines import sylph as sylph_mod

                sylph_db = Path(args.refs_dir) / "sylph"
                if sylph_mod.available(sylph_db):
                    with reporter.stage("whole-genome profile"):
                        if taxonomy is not None and taxonomy.host is not None:
                            g1, g2 = taxonomy.host.nonhost_r1, taxonomy.host.nonhost_r2
                        else:
                            g1 = sample_input.mates[0].path
                            g2 = sample_input.mates[1].path if len(sample_input.mates) > 1 else None
                        genome_profile = sylph_mod.run_sylph(
                            sample=str(sample_input.sample), r1=g1, r2=g2,
                            db_dir=sylph_db, work_dir=out_dir / "genome",
                            threads=args.threads,
                        )
                        if genome_profile is not None:
                            reporter.record(
                                f"    {len(genome_profile.confident_hits)} species clusters "
                                f"against {sylph_mod.CATALOGUE_SIZE:,} in GTDB"
                                + (" (cached)" if genome_profile.cached else "")
                            )
                else:
                    reporter.record("    whole-genome profile: tool or database not installed; lane not assessed")
            except Exception as exc:  # noqa: BLE001 - an optional lane never loses a run
                reporter.warn(f"whole-genome profile failed: {exc}")
                genome_profile = None
            timings["whole-genome profile"] = time.monotonic() - t1

        # Pathogen branch (spec v5.0). Runs on host-filtered reads when the
        # host step ran, and on the raw pair otherwise: a missing host filter
        # costs specificity that the branch's own decoy masking recovers, but
        # a missing branch costs the reader the entire question.
        if not args.no_pathogens:
            t1 = time.monotonic()
            try:
                with reporter.stage("pathogen screen"):
                    if taxonomy is not None and taxonomy.host is not None:
                        p1, p2 = taxonomy.host.nonhost_r1, taxonomy.host.nonhost_r2
                    else:
                        p1 = sample_input.mates[0].path
                        p2 = (
                            sample_input.mates[1].path
                            if len(sample_input.mates) > 1
                            else None
                        )
                    pathogens = run_pathogen_branch(
                        sample_id=sample_input.sample,
                        r1=p1,
                        r2=p2,
                        catalog=load_pathogen_catalog(Path("pathogens/catalog")),
                        bundle_dir=Path(args.refs_dir) / "pathogens" / "bundle",
                        determinant_dir=Path(args.refs_dir) / "pathogens" / "determinants",
                        host_fasta=Path(args.refs_dir) / "host" / "grch38.fna.gz",
                        out_dir=out_dir,
                        # Was hardcoded True and then reported back per target
                        # as `qc_pass`, so the report carried a QC verdict that
                        # had never been evaluated.
                        qc_passed=(qc is None or not qc.failed),
                        threads=args.threads,
                        progress=reporter.step,
                        strict_cdiff=bool(getattr(args, "strict_cdiff", False)),
                    )
                    reporter.info(
                        f"{pathogens.pathogen_count} supported pathogen finding(s); "
                        f"{pathogens.coverage.assessed:,} targets assessed, "
                        f"{pathogens.coverage.not_assessed:,} not assessed"
                    )
            except DependencyError as exc:
                reporter.warn(f"pathogen screen skipped: {exc}")
            except OpenBiotaError as exc:
                reporter.warn(f"pathogen screen failed: {exc}")
            timings["pathogen screen"] = time.monotonic() - t1

        # Estimated stool-microbiome chronological age (research use), spec §9.
        # Runs on the MetaPhlAn 3 species profile only — the frozen bundle's
        # namespace — and fails closed without a known adult age in participant mode.
        if taxonomy is not None and not args.no_age:
            t1 = time.monotonic()
            with reporter.stage("estimated chronological age (research use)"):
                bundle = load_age_bundle(Path(args.refs_dir))
                if bundle is None:
                    reporter.warn(
                        f"no age model bundle at {age_bundle_dir(Path(args.refs_dir))} — run "
                        "`openbiota age-model train` to build it; the age page is omitted"
                    )
                else:
                    age_meta = dict(bundle.metadata)
                    age_result = predict_age(
                        bundle, taxonomy.species,
                        actual_age=args.subject_age, age_known=args.subject_age is not None,
                        mode=args.mode, profiler_index=taxonomy.index,
                        # MetaPhlAn 3 does not write a read count into its profile;
                        # the depth gate then uses the non-host pair count (x2 mates).
                        n_reads=(
                            taxonomy.n_reads_processed
                            or ((taxonomy.host.total_pairs - taxonomy.host.host_pairs) * len(sample_input.mates)
                                if taxonomy.host is not None and taxonomy.host.total_pairs else None)
                            or (known_total_pairs * len(sample_input.mates) if known_total_pairs else None)
                        ),
                        unknown_percent=taxonomy.unknown_percent,
                        country=args.subject_country,
                        qc_failed=False,
                    )
                    reporter.record(
                        "    age: " + (
                            f"{age_result.status_line}, estimate {age_result.predicted_chronological_age_years} y "
                            f"(95% {age_result.prediction_interval_95[0]}–{age_result.prediction_interval_95[1]})"
                            if age_result.status == "scored"
                            else f"{age_result.status_line} ({age_result.reason})"
                        )
                    )
            timings["chronological age"] = time.monotonic() - t1

        t0 = time.monotonic()
        with reporter.stage("profile similarity"):
            similarity = run_similarity_stage(
                profile_set=load_profile_set(args.profiles_dir),
                profile_names=_split_list(args.profile_names),
                tally=result,
                taxonomy=taxonomy,
                cohort_path=Path(args.taxonomic_cohort),
                functional_samples_path=Path(args.functional_cohort),
                subject=subject,
                match_on=[m.strip() for m in args.match.split(",") if m.strip()],
                reporter=reporter,
                include_opt_in=args.include_opt_in,
                profile_validation=load_optional_json(Path(args.profile_validation)),
                taxon_groups=load_taxon_group_set(Path(args.taxa_dir)),
                # The skin panels have their own depth floor, stated in pairs.
                usable_read_pairs=known_total_pairs,
                qc_passed=qc is None or not qc.failed,
                mode=args.mode,
            )
        timings["profile similarity"] = time.monotonic() - t0

    # The age model sees species composition only, and on its own regresses
    # towards its adult training mean. GMWI2 and the scored disease patterns
    # are two health readings it never sees, and both carry age information;
    # applying them here — after similarity has run — moves the estimate
    # inside the model's own 90% interval. Order matters: this needs both
    # stages to have completed, and it must not run twice.
    if age_result is not None and similarity is not None:
        shape = similarity.shape
        age_result = adjust_for_community_health(
            age_result,
            gmwi2=None if shape is None else shape.dysbiosis_score,
            # Only patterns that actually scored: an abstention is not a
            # favourable reading and must not take a year off.
            profile_percentiles=[
                r.combined_percentile
                for r in similarity.results
                if r.reportable and r.combined_percentile is not None
            ],
        )
        if age_result.health_adjustment is not None:
            d = age_result.health_adjustment
            reporter.info(
                f"biological age {age_result.health_adjusted_age_years} y "
                f"(model {d['model_point_years']}, "
                f"{d['disease_years']:+g} y from disease patterns, "
                f"{d['gmwi2_years']:+g} y from GMWI2"
                + (f", held at the {d['clamped_to']}" if d["clamped"] else "")
                + ")"
            )

    # ---- report ----------------------------------------------------------- #
    selected_results = [
        r for r in result.panels if r.panel.name in {p.name for p in selected_panels}
    ]
    diagnostics = run_diagnostics(tally=result, qc=qc, selected=selected_results)

    elapsed = time.monotonic() - started
    run_meta: dict[str, Any] = {
        "openbiota version": __version__,
        "started": started_at,
        "wall clock": human_duration(elapsed),
        "diamond": database.diamond_version.splitlines()[0],
        "diamond flags": _flags_for_display(searches[0].command),
        "sensitivity": args.sensitivity,
        "threads": args.threads,
        "block size": args.block_size,
        "index chunks": args.index_chunks,
        "e-value": args.evalue,
        "reference fingerprint": database.fingerprint,
        "reference proteins": f"{database.n_sequences:,} ({database.total_residues():,} aa)",
        "reference built": database.built_at,
        "panels reported": ", ".join(p.name for p in selected_panels),
        "taxonomic engine": (
            f"{similarity.taxonomy.family} {similarity.taxonomy.profiler_version}, "
            f"database {similarity.taxonomy.index}"
            if similarity is not None and similarity.taxonomy is not None
            else "not run"
        ),
        "taxonomic reference": (
            f"cMD {similarity.cohort_manifest.get('snapshot')} "
            f"n={similarity.cohort_manifest.get('n_samples'):,} "
            f"manifest {similarity.cohort_manifest.get('manifest_id')}"
            if similarity is not None and similarity.cohort_manifest
            else "none"
        ),
        "profiles scored": (
            ", ".join(f"{r.profile.name} v{r.profile.version}" for r in similarity.results)
            if similarity is not None
            else "none"
        ),
        "subsample": args.subsample or "no (full input)",
        "output directory": str(out_dir),
        "stage timings": ", ".join(f"{k} {human_duration(v)}" for k, v in timings.items()),
    }

    summary_plain = render_summary(
        sample=sample_input.sample,
        tally=result,
        selected=selected_results,
        qc=qc,
        profile=profile,
        diagnostics=diagnostics,
        run_meta=run_meta,
        style=Style(False),
        similarity=similarity,
    )
    results_json = build_results_json(
        sample=sample_input.sample,
        tally=result,
        selected=selected_results,
        qc=qc,
        profile=profile,
        diagnostics=diagnostics,
        run_meta={k: str(v) for k, v in run_meta.items()},
        reference_entries={k: v.to_json() for k, v in sorted(database.entries.items())},
    )
    # ---- cohort comparison + PDF ------------------------------------------ #
    ranges: ReferenceRanges | None = None
    if Path(args.reference_ranges).is_file():
        try:
            ranges = ReferenceRanges.load(Path(args.reference_ranges))
        except (OSError, ValueError, KeyError) as exc:
            reporter.warn(f"reference ranges unusable ({exc}); percentiles unavailable")
    else:
        reporter.warn(
            f"no reference ranges at {args.reference_ranges} — run `openbiota cohort` to build "
            "them, otherwise results have no population context"
        )

    # The organisms behind each reading come from the fragment tables joined
    # to the sample's organism list. That list is assembled in full later
    # (with strains); here a strain-less copy from what is already in hand
    # is enough to say whether a driver is present and what class it is.
    _driver_inv = None
    try:
        from openbiota import inventory as _inv_mod

        _driver_inv = _inv_mod.from_results(
            {
                "extended_catalogue": (
                    extended.to_json() if extended is not None else {}
                ),
                "genome_profile": (
                    genome_profile.to_json() if genome_profile is not None else {}
                ),
                "detection": detection_block,
            },
            scoring_species=(
                similarity.taxonomy.species
                if similarity is not None and similarity.taxonomy is not None
                else None
            ),
            ranked_rows=(
                similarity.community.species
                if similarity is not None and getattr(similarity, "community", None) is not None
                else ()
            ),
        )
    except Exception as exc:  # noqa: BLE001 - drivers are an enrichment, never a blocker
        reporter.warn(f"driver inventory unavailable: {exc}")

    report_rows = build_rows(
        selected_results, ranges, sample_dir=out_dir, inventory=_driver_inv,
    )
    results_json["reference_comparison"] = reference_json(report_rows, ranges)
    # Every displayed value, written once, so the report is drawn from the
    # file rather than from a second derivation of the same panels. A
    # renderer that recomputes can disagree with the JSON it was handed, and
    # the web interface that will read this file later has to see exactly
    # what the PDF saw.
    results_json["report_rows"] = [row.to_json() for row in report_rows]
    results_json["metabolite_drivers"] = {
        r.panel: (r.drivers.to_json() if r.drivers is not None else None) for r in report_rows
    }
    if similarity is not None:
        results_json["profile_similarity"] = similarity.to_json()
    if age_result is not None:
        results_json["microbiome_age"] = age_result.to_json()
    if pathogens is not None:
        results_json["pathogens"] = pathogens.to_json()
    results_json["subject_context"] = context_ledger.to_json()

    # ---- sequencing quality gates (spec 4.2) ------------------------------ #
    host = similarity.taxonomy.host if similarity is not None and similarity.taxonomy else None
    classified = (
        None if similarity is None or similarity.taxonomy is None
        else max(0.0, 1.0 - similarity.taxonomy.unknown_percent / 100.0)
    )
    sampled_mate = qc.mates[0] if qc is not None and qc.mates else None
    gates = assemble_gates(
        input_pairs=qc.read_pairs if qc is not None else (fastp_result.pairs_in if fastp_result else None),
        fastp=fastp_result,
        host_fraction=None if host is None else host.host_fraction,
        host_status=None if host is None else host.status,
        nonhost_pairs=None if host is None else host.nonhost_pairs,
        classified_fraction=classified,
        rpob_fragments=result.normalizer.fragments,
        sampled_duplicate_fraction=None if sampled_mate is None else sampled_mate.duplicate_fraction,
        sampled_gc_percent=None if sampled_mate is None else sampled_mate.gc_percent,
        sampled_mean_length=None if sampled_mate is None else sampled_mate.mean_length,
        provenance=provenance,
    )
    results_json["sequencing_quality"] = {
        "gates": gates.to_json(),
        "fastp": None if fastp_result is None else fastp_result.to_json(),
    }
    if extended is not None:
        results_json["extended_catalogue"] = extended.to_json(
            mpa3_species=(
                similarity.taxonomy.species
                if similarity is not None and similarity.taxonomy is not None
                else None
            )
        )
    results_json["detection"] = detection_block
    if genome_profile is not None:
        results_json["genome_profile"] = genome_profile.to_json()
    else:
        results_json["genome_profile"] = {
            "status": "not_assessed",
            "what_this_is": (
                "Whole-genome profiling against GTDB did not run for this sample: the "
                "tool or its database was not installed. Nothing here says which "
                "organisms that lane would have found."
            ),
        }
    # "P. copri" is 13 species-level clades, 13-21% divergent from each other,
    # and the extended lane already resolves every one of them. Reporting the
    # composition turns the rheumatoid-arthritis profile's central ambiguity
    # into a measured statement — usually a definite absence across the whole
    # complex, which a species-level number could never assert.
    from openbiota import copri

    results_json["copri_complex"] = copri.resolve(
        (results_json.get("extended_catalogue") or {}).get("sgbs")
        if extended is not None
        else None
    )

    # Strain resolution: dominant population fingerprints per organism, from
    # the marker lane. Present only when that lane has been run for this
    # sample; absent means "not yet run", never "no strains" (spec §3.2).
    from openbiota.resolution import strains as strain_mod

    # These last stages read artefacts produced by separate passes - a
    # bz2-compressed marker file, a locus index. A truncated one raises
    # inside json/bz2, and because the reporter buffers its lines until the
    # very end, an uncaught error here lost the summary, results.json and
    # run.log for an analysis that had otherwise completed. Each is now
    # wrapped: the block records why it is missing and the run still writes
    # everything else.
    try:
        results_json["strain_resolution"] = _strain_block(
            out_dir,
            str(sample_input.sample),
            results_json.get("extended_catalogue"),
            strain_mod,
        )
    except Exception as exc:  # noqa: BLE001 - never lose a run for this block
        reporter.warn(f"strain resolution unavailable: {exc}")
        results_json["strain_resolution"] = {
            "status": "stage_error",
            "error": str(exc),
            "what_this_is": (
                "Strain resolution did not complete for this sample. That is "
                "different from finding no strains: nothing was resolved, so "
                "nothing here says which populations are present."
            ),
        }

    # One organism inventory, merged from every profiler that ran, under
    # current taxonomic names. This is what the report's organism listings
    # read: reporting a single catalogue discards real detections, and
    # reporting two side by side makes the reader reconcile them by hand.
    # Percentiles stay bound to the catalogue the reference cohort was
    # processed on - see openbiota/inventory.py for why that line is hard.
    from openbiota import inventory as inventory_mod

    try:
        _inv = inventory_mod.from_results(
            results_json,
            scoring_species=(
                similarity.taxonomy.species
                if similarity is not None and similarity.taxonomy is not None
                else None
            ),
            ranked_rows=(
                similarity.community.species
                if similarity is not None
                and getattr(similarity, "community", None) is not None
                else ()
            ),
        )
        results_json["organism_inventory"] = _inv.to_json()
        # Competitive confirmation (spec 0.8.4 §4C): every organism the
        # installed baseline never saw, and every single-method call, has to
        # win its reads against its detected relatives before it is more
        # than provisional. Rejected calls stay in the file as rejected.
        t_conf = time.monotonic()
        try:
            from openbiota.expansion import confirm as _confirm

            if _confirm.available() and _confirm.select_candidates(_inv):
                with reporter.stage("competitive confirmation"):
                    if taxonomy is not None and taxonomy.host is not None:
                        c1, c2 = taxonomy.host.nonhost_r1, taxonomy.host.nonhost_r2
                    else:
                        c1 = sample_input.mates[0].path
                        c2 = sample_input.mates[1].path if len(sample_input.mates) > 1 else None
                    confirmation = _confirm.run_confirmation(
                        sample=str(sample_input.sample), inventory=_inv, r1=c1, r2=c2,
                        work_dir=out_dir / "confirmation", threads=args.threads,
                    )
                    results_json["organism_inventory"] = _confirm.apply_verdicts(
                        results_json["organism_inventory"], confirmation)
                    results_json.setdefault("detection", {})["confirmation"] = confirmation
                    vs_ = confirmation.get("verdicts") or []
                    reporter.record(
                        f"    {confirmation.get('n_tested', 0)} candidates tested against "
                        f"{confirmation.get('n_genomes_in_competition', 0)} genomes: "
                        f"{sum(1 for v in vs_ if v['status'] == 'supported')} supported, "
                        f"{sum(1 for v in vs_ if v['status'] == 'provisional')} provisional, "
                        f"{sum(1 for v in vs_ if v['status'] == 'not_detected')} rejected, "
                        f"{sum(1 for v in vs_ if v['status'] == 'pending')} pending a genome"
                    )
                    # The in-memory inventory must agree with the file.
                    _inv = inventory_mod.from_json(results_json["organism_inventory"]) or _inv
            else:
                results_json.setdefault("detection", {})["confirmation"] = {
                    "status": "completed", "n_candidates": 0, "verdicts": [],
                    "note": "no candidate needed confirmation" if _confirm.available() else "bowtie2/samtools not installed"}
                if isinstance(results_json.get("organism_inventory"), dict):
                    inventory_mod.unify_shares(results_json["organism_inventory"], None)
                    _inv = inventory_mod.from_json(results_json["organism_inventory"]) or _inv
                    results_json["organism_inventory"] = _inv.to_json()
        except Exception as exc:  # noqa: BLE001 - never lose a run for this stage
            reporter.warn(f"competitive confirmation failed: {type(exc).__name__}: {exc}")
            results_json.setdefault("detection", {})["confirmation"] = {"status": "failed", "reason": f"{type(exc).__name__}: {exc}"}
        timings["competitive confirmation"] = time.monotonic() - t_conf

        # Comparative strain analysis (spec 0.8.4 §5): StrainPhlAn on the
        # Jan26 SAM against member genomes of each organism's R232 species.
        # The Jun23 consensus block above is the preserved baseline; this is
        # the comparative workflow it never was.
        t_sp = time.monotonic()
        try:
            from openbiota.expansion import strainphlan as _sp

            jan26_sam = out_dir / "strain" / f"{sample_input.sample}.jan26.markers.sam.bz2"
            if _sp.available() and jan26_sam.is_file():
                with reporter.stage("comparative strain analysis (StrainPhlAn 4 / Jan26)"):
                    results_json["strain_analysis"] = _sp.run(
                        sample=str(sample_input.sample), sam=jan26_sam, db_dir=Path(args.refs_dir) / "metaphlan4_db",
                        work_dir=out_dir / "strain", threads=args.threads, inventory=_inv,
                    )
                    sa = results_json["strain_analysis"]
                    reporter.record(
                        f"    {sa.get('n_compared', 0)} organisms placed against reference genomes, "
                        f"{sa.get('n_unresolved', 0)} unresolved, {len(sa.get('deferred') or [])} deferred to the next run"
                    )
                    if isinstance(results_json.get("organism_inventory"), dict) and \
                            _sp.attach_to_inventory(results_json["organism_inventory"], sa):
                        _inv = inventory_mod.from_json(results_json["organism_inventory"]) or _inv
            else:
                results_json["strain_analysis"] = {
                    "status": "not_assessed",
                    "reason": ("StrainPhlAn externals missing" if not _sp.available() else "no Jan26 marker SAM for this sample"),
                }
        except Exception as exc:  # noqa: BLE001
            reporter.warn(f"comparative strain analysis failed: {type(exc).__name__}: {exc}")
            results_json["strain_analysis"] = {"status": "failed", "reason": f"{type(exc).__name__}: {exc}"}
        timings["strain analysis"] = time.monotonic() - t_sp

        # BacDive enrichment (spec 0.8.4 §6): culture-collection identities
        # and reference-strain traits for the named species, most abundant
        # first, inside a time budget. Cached, so re-runs cost nothing.
        t_bd = time.monotonic()
        try:
            from openbiota.expansion import bacdive as _bacdive

            with reporter.stage("BacDive enrichment (API v2)"):
                named = [o.species for o in _inv.organisms if not o.unnamed]
                results_json["bacdive"] = _bacdive.enrich(named, budget_s=240.0)
                reporter.record(
                    f"    {results_json['bacdive']['n_found']} of {results_json['bacdive']['n_species_queried']} "
                    f"named species have BacDive records"
                )
        except Exception as exc:  # noqa: BLE001
            reporter.warn(f"BacDive enrichment failed: {type(exc).__name__}: {exc}")
            results_json["bacdive"] = {"status": "failed", "reason": f"{type(exc).__name__}: {exc}"}
        timings["bacdive"] = time.monotonic() - t_bd

        # Targeted assembly (spec 0.8.4 §4D): runs when SingleM left lineages
        # unresolved at reconstructable coverage. A data condition, not a flag.
        t_as = time.monotonic()
        try:
            from openbiota.expansion import assembly as _assembly

            trig = ((((results_json.get("detection") or {}).get("lanes") or {}).get("singlem_globdb") or {})
                    .get("summary") or {}).get("assembly_triggers") or []
            if trig and not getattr(args, "no_assembly", False):
                with reporter.stage(f"targeted assembly ({len(trig)} unresolved lineages; metaSPAdes)"):
                    if taxonomy is not None and taxonomy.host is not None:
                        a1, a2 = taxonomy.host.nonhost_r1, taxonomy.host.nonhost_r2
                    else:
                        a1 = sample_input.mates[0].path
                        a2 = sample_input.mates[1].path if len(sample_input.mates) > 1 else None
                    refs_fetched = sorted(Path("refs/genomes").glob("*.fna.gz"))[:400]
                    # Targeted: only the reads Kraken/UHGG could not place at
                    # species level are assembled (assembly.select_unplaced_reads).
                    k_assign = out_dir / "kraken" / f"{sample_input.sample}.kraken2.assignments.tsv.gz"
                    k_db = Path(args.refs_dir) / "kraken2" / "uhgg_v2.0.2"
                    results_json["assembly"] = _assembly.run(
                        sample=str(sample_input.sample), r1=a1, r2=a2, work_dir=out_dir / "assembly",
                        threads=args.threads, triggers=trig, reference_genomes=refs_fetched,
                        kraken_assignments=k_assign if k_assign.is_file() else None,
                        kraken_db=k_db if (k_db / "taxonomy" / "nodes.dmp").is_file() else None,
                    )
                    asm = results_json["assembly"]
                    if asm.get("status") == "completed":
                        st_ = asm["assembly"]
                        sel_ = asm.get("reads") or {}
                        if sel_.get("n_selected_pairs") is not None:
                            reporter.record(f"    assembled {sel_['n_selected_pairs']:,} read pairs "
                                            f"({sel_.get('selected_fraction', 0):.1%}) that no gut catalogue placed at species level")
                        reporter.record(f"    {st_['n_contigs']:,} contigs >= {asm['min_contig']} bp, N50 {st_['n50']:,}, "
                                        f"{sum(1 for t in asm['triggers'] if t['n_marker_families_on_contigs'])} lineages with marker-bearing contigs")
            else:
                results_json["assembly"] = {"status": "not_triggered" if not trig else "skipped",
                                            "reason": "no unresolved lineage reached the assembly coverage" if not trig else "--no-assembly"}
        except Exception as exc:  # noqa: BLE001
            reporter.warn(f"targeted assembly failed: {type(exc).__name__}: {exc}")
            results_json["assembly"] = {"status": "failed", "reason": f"{type(exc).__name__}: {exc}"}
        timings["assembly"] = time.monotonic() - t_as
        # Section 4 reads the inventory for what each group contains. The
        # percentile stays on the reference catalogue's namespace; detection
        # and composition of every member come from every lane.
        if similarity is not None and getattr(similarity, "community", None) is not None:
            from openbiota.community import enrich_with_inventory

            enrich_with_inventory(similarity.community, _inv)
            results_json["profile_similarity"] = similarity.to_json()
            n_exp = sum(g.n_detected_expansion_only for g in similarity.community.groups)
            if n_exp:
                reporter.record(f"    {n_exp} group member(s) detected only by the expanded lanes")
        # The pathogen screen aligns to its own bundle, so a relative it
        # carries no reference for has its shared sequence credited to the
        # nearest target. Reconcile every bacterial call with the inventory,
        # which mapped the same reads competitively against the relatives it
        # actually found, so the two sections cannot contradict each other.
        if isinstance(results_json.get("pathogens"), dict):
            try:
                from openbiota.pathogens import reconcile as _reconcile

                _rec = _reconcile.reconcile(
                    results_json["pathogens"], _inv,
                    read_length_bp=((results_json.get("sequencing_quality") or {}).get("fastp") or {}).get("read1_mean_length"),
                    rejected=(results_json.get("organism_inventory") or {}).get("rejected") or [],
                )
                if _rec.get("n_shared_sequence"):
                    reporter.record(f"    pathogen screen: {_rec['n_shared_sequence']} bacterial call(s) are shared "
                                    f"sequence from relatives the inventory found; {_rec['n_found']} agree with it")
            except Exception as exc:  # noqa: BLE001 - a reconciliation pass, never fatal
                reporter.warn(f"pathogen reconciliation failed: {type(exc).__name__}: {exc}")
        # One verdict per organism: class, description basis, and whether
        # its level is a concern. Computed here so the report and any
        # downstream consumer read the same judgement.
        from openbiota import organisms as organisms_mod

        _vs = organisms_mod.verdicts(_inv.organisms)
        results_json["organism_verdicts"] = {
            "composition_percent_of_reads": organisms_mod.composition(_vs),
            "n_flagged": sum(1 for v in _vs if v.flagged),
            "n_issues": sum(1 for v in _vs if v.is_issue),
            "verdicts": [v.to_json() for v in _vs],
            "vocabulary": {
                c: organisms_mod.CLASS_MEANING[c] for c in organisms_mod.CLASSES
            },
        }
    except Exception as exc:  # noqa: BLE001 - never lose a run for this block
        reporter.warn(f"organism inventory unavailable: {exc}")
        results_json["organism_inventory"] = {"status": "stage_error", "error": str(exc)}

    # Targeted locus mapping: subtype discrimination and locus architecture at
    # nucleotide level. This is the stage that can contradict a translated
    # protein hit -- a gene-family match at high amino-acid identity can be to
    # a homolog whose gene shares no nucleotide similarity with the real
    # target, and only a nucleotide check can tell.
    from openbiota.resolution import loci as loci_mod

    try:
        locus_calls = loci_mod.run_all(
            sample_id=str(sample_input.sample),
            reads=[p for p in (r1, r2) if p is not None],
            threads=args.threads,
        )
        results_json["locus_resolution"] = {
            "targets": [c.to_json() for c in locus_calls],
            "what_this_is": (
                "Reads mapped onto exact accession-versioned nucleotide references for two "
                "questions a protein search cannot answer: which allele of a toxin is present, "
                "and whether a multi-gene biosynthetic island is present as a locus rather than "
                "as scattered fragments."
            ),
        }
    except Exception as exc:  # noqa: BLE001 - never lose a run for this block
        reporter.warn(f"locus resolution unavailable: {exc}")
        results_json["locus_resolution"] = {
            "targets": [],
            "status": "stage_error",
            "error": str(exc),
            "what_this_is": (
                "Nucleotide locus mapping did not complete. No allele or locus "
                "call was made, which is different from making one and finding "
                "nothing."
            ),
        }

    # The completeness census: every registered resolution target with the
    # state it actually reached. A partial run cannot present itself as a
    # complete screen (spec §12.2, acceptance 1-2).
    try:
        results_json["resolution_census"] = _resolution_census(
            str(sample_input.sample), results_json
        )
    except Exception as exc:  # noqa: BLE001 - never lose a run for this block
        reporter.warn(f"resolution census unavailable: {exc}")
        results_json["resolution_census"] = {
            "status": "stage_error",
            "error": str(exc),
            "what_this_is": (
                "The completeness census did not build, so this run cannot "
                "show which registered targets were reached."
            ),
        }

    # Biofilm-related potential (spec v08.0 §9). Two independent axes, never
    # combined. The engine needs the MetaPhlAn 3 profile, which shares a
    # taxonomy namespace with the reference cohort; if it is missing the block
    # still renders, reporting the incompatibility rather than ranking a
    # MetaPhlAn 4 abundance against a MetaPhlAn 3 reference.
    try:
        from openbiota.biofilm import engine as _bf_engine

        results_json["biofilm"] = _bf_engine.analyze(out_dir, results_json)
    except Exception as exc:  # noqa: BLE001 - never fail a run for this block
        reporter.warn(f"biofilm module unavailable: {exc}")
        results_json["biofilm"] = {
            "schema_version": "openbiota.biofilm/1.1",
            "sample_id": str(sample_input.sample),
            "status": "module_error",
            "error": str(exc),
            "cards": [],
            "limits": [
                "The biofilm module did not complete. No biofilm result is "
                "reported for this sample, which is different from a result of "
                "no biofilm-related findings."
            ],
        }

    # The mycobiome module (spec v08.2): EukDetect2 marker lane, competitive
    # whole-genome fungal lane over every read pair, strain attempt for every
    # supported fungus, evidence registry and the experimental MHS-E1 score.
    # Runs when the fungal reference lock exists; otherwise not assessed.
    try:
        from openbiota.mycobiome import engine as _myco_engine
        from openbiota.mycobiome import references as _myco_refs

        if _myco_refs.available(Path(args.refs_dir) / "mycobiome"):
            with reporter.stage("mycobiome"):
                if taxonomy is not None and taxonomy.host is not None:
                    m1, m2 = taxonomy.host.nonhost_r1, taxonomy.host.nonhost_r2
                else:
                    m1 = sample_input.mates[0].path
                    m2 = sample_input.mates[1].path if len(sample_input.mates) > 1 else None
                sketch = out_dir / "genome" / "sketches" / f"{sample_input.sample}.paired.sylsp"
                results_json["mycobiome"] = _myco_engine.analyze(
                    sample=str(sample_input.sample), r1=m1, r2=m2, results=results_json, out_dir=out_dir,
                    refs_dir=Path(args.refs_dir), threads=args.threads,
                    read_sketch=sketch if sketch.is_file() else None, log=reporter.record,
                )
        else:
            results_json["mycobiome"] = {
                "schema_version": "openbiota.mycobiome.v3", "sample_id": str(sample_input.sample),
                "analysis_status": "not_assessed",
                "limits": ["Fungal reference lock not built; run `make mycobiome-refs`."],
            }
    except Exception as exc:  # noqa: BLE001 - never fail a run for this module
        reporter.warn(f"mycobiome module unavailable: {exc}")
        results_json["mycobiome"] = {
            "schema_version": "openbiota.mycobiome.v3", "sample_id": str(sample_input.sample),
            "analysis_status": "failed", "error": str(exc),
            "limits": ["The mycobiome module did not complete. No fungal result is reported, which is "
                       "different from a result of no fungi."],
        }


    reporter.record(
        f"    QC gates: {gates.overall}"
        + (f", {gates.usable_pairs:,} usable non-host pairs" if gates.usable_pairs is not None else "")
    )
    if gates.depth_adequate is False:
        reporter.warn(
            f"below the prespecified {MIN_USABLE_NONHOST_PAIRS:,} usable non-host pairs: "
            "taxonomic non-detections are reported as missing, not absent"
        )

    # ---- organism findings + intervention evidence (spec v04 §5–8) -------- #
    plan: ActionPlan | None = None
    t0 = time.monotonic()
    # Share of each panel's gene-family fragments that a scope rule demoted to
    # background: how much of the family sits in organisms that carry the fold
    # without the pathway. Reported beside the capacity figure.
    out_of_scope: dict[str, float] = {}
    for pr in result.panels:
        oos = sum(d.fragments for d in pr.decoys if d.entry_id.endswith(OUT_OF_SCOPE_SUFFIX))
        if oos and (pr.accepted_fragments + oos):
            out_of_scope[pr.panel.name] = oos / (pr.accepted_fragments + oos)
    with reporter.stage("organism findings and intervention evidence"):
        try:
            plan = build_action_plan(
                sample=sample_input.sample, rows=report_rows, similarity=similarity, gates=gates,
                subject_age=args.subject_age, subject_sex=args.subject_sex,
                subject_country=args.subject_country,
                medications=_split_list(args.medications) if args.medications else None,
                mode=args.mode, taxa_dir=Path(args.taxa_dir), reporter=reporter,
                out_of_scope_fractions=out_of_scope, context_ledger=context_ledger,
            )
        except OpenBiotaError as exc:
            reporter.warn(f"intervention evidence unavailable: {exc}")
        except Exception as exc:  # noqa: BLE001 — a report layer, never fatal to the run
            reporter.warn(
                f"organism findings / evidence failed: {type(exc).__name__}: {exc}\n"
                + traceback.format_exc(limit=6)
            )
    timings["findings and evidence"] = time.monotonic() - t0
    if plan is not None:
        results_json["findings_and_evidence"] = plan.to_json()

    # ---- v0.8.3 extension (spec 0.8.3 §3.2) ------------------------------ #
    # Runs after the evidence plan is attached: the A12 planner consolidates
    # those cards, and the A14 register records what the other views assumed.
    # Additive by construction: every new measurement is written under one new
    # key and no protected object is touched. A failure here costs the new
    # views and nothing else, which is why it cannot raise into the run.
    t0 = time.monotonic()
    with reporter.stage("report extension"):
        try:
            from openbiota.extension import engine as _ext_engine

            _extension = _ext_engine.analyze(
                results_json, mode="full-extension" if args.simulate != "never" else "cached",
                sample=str(sample_input.sample),
                simulate=getattr(args, "simulate", "auto"),
                lp_workers=getattr(args, "lp_workers", None),
                progress=reporter.record,
            )
            _ext_engine.attach(results_json, _extension)
            _ext_files = _ext_engine.write_outputs(_extension, out_dir)
            reporter.record(
                f"    {len(_extension.metrics)} additional measurements, "
                f"{len(_extension.unavailable)} capabilities unavailable "
                f"({len(_extension.manifest.get('engineering_debt') or [])} not yet implemented)"
            )
            reporter.record("    wrote " + ", ".join(sorted(_ext_files)))
        except Exception as exc:  # noqa: BLE001 - the extension never fails a run
            reporter.warn(f"report extension unavailable: {type(exc).__name__}: {exc}")
    timings["report extension"] = time.monotonic() - t0

    # The last word on the run: did every stage produce its result? Judged
    # from the results object itself, recorded in it, and printed at the end.
    # A stage that warned and moved on shows up here as MISSING, and the run
    # exits non-zero for it - the report must never quietly go out with a
    # section that says it did not run.
    from openbiota.completeness import audit as _audit_completeness  # noqa: PLC0415

    completeness = _audit_completeness(results_json)
    results_json.setdefault("run", {})["completeness"] = completeness.to_json()
    # Spec 0.8.4 §7 reproducibility gate: actual peak RAM, disk and runtime
    # are recorded, not estimated. ru_maxrss is bytes on macOS, KB on Linux.
    try:
        import resource as _resource
        import sys as _sys

        _scale = 1 if _sys.platform == "darwin" else 1024
        _self = _resource.getrusage(_resource.RUSAGE_SELF).ru_maxrss * _scale
        _kids = _resource.getrusage(_resource.RUSAGE_CHILDREN).ru_maxrss * _scale
        _disk = sum(p.stat().st_size for p in out_dir.rglob("*") if p.is_file())
        results_json["run"]["resources"] = {
            "peak_rss_gb_pipeline": round(_self / 1e9, 2),
            "peak_rss_gb_largest_child": round(_kids / 1e9, 2),
            "results_dir_gb": round(_disk / 1e9, 2),
            "wall_s": round(time.monotonic() - started, 1),
            "stage_timings_s": {k: round(v, 1) for k, v in timings.items()},
            "threads": args.threads,
            "note": "peak resident set of the pipeline process and of its largest single child (the sum over "
                    "children is not observable after the fact); disk is the sample's results directory",
        }
    except Exception:  # noqa: BLE001 - accounting never fails a run
        pass

    summary_path, json_path = write_outputs(
        out_dir=out_dir, summary_text=summary_plain, results=results_json
    )

    pdf_path: Path | None = None
    pdf_error: str | None = None
    if not args.no_pdf:
        t0 = time.monotonic()
        with reporter.stage("PDF report"):
            try:
                pdf_path = build_pdf(
                    path=out_dir / f"{sample_input.sample}_report.pdf",
                    results=results_json,
                    rows=report_rows,
                    profile=None if profile is None else profile.to_json(),
                    validation=load_optional_json(Path(args.validation)),
                    cohort=None if ranges is None else ranges.to_json(),
                    cohort_note=cohort_note(ranges),
                    depth=load_optional_json(out_dir / "depth_check.json"),
                    similarity=similarity,
                    profile_validation=load_optional_json(Path(args.profile_validation)),
                    plan=plan,
                    age=age_result,
                    age_meta=age_meta,
                    gates=gates,
                    mode=args.mode,
                    subject_age=args.subject_age,
                )
            except ImportError:
                reporter.warn("reportlab is not installed; skipping the PDF report")
            except Exception as exc:  # noqa: BLE001 — keep summary/JSON, but fail loudly
                reporter.error(
                    f"PDF report failed: {type(exc).__name__}: {exc}\n"
                    + traceback.format_exc(limit=6)
                )
                pdf_error = f"{type(exc).__name__}: {exc}"
        timings["pdf report"] = time.monotonic() - t0

    # The whole run, start to finish, PDF included. `elapsed` above was taken
    # before the PDF stage and is the analysis time; this is what the reader
    # of the log wants to know. It goes in the header and, in every mode, as
    # the log's last line — emitted *before* the buffer is flushed to disk,
    # which is the only way it can reach the file at all.
    total = time.monotonic() - started
    produced = [summary_path.name, json_path.name, "run.log"]
    if pdf_path is not None:
        produced.insert(0, pdf_path.name)
    results_json["run"]["total wall clock"] = human_duration(total)
    results_json["run"]["total_s"] = round(total, 1)
    json_path.write_text(json.dumps(results_json, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    if pdf_error is not None:
        reporter.error(
            f"finished in {human_duration(total)} WITHOUT a PDF ({pdf_error}) — "
            f"{', '.join(produced)} in {out_dir}"
        )
    elif not completeness.complete:
        from openbiota.completeness import describe as _describe_completeness  # noqa: PLC0415

        reporter.error(
            f"INCOMPLETE RUN in {human_duration(total)}: {len(completeness.missing)} of "
            f"{len(completeness.stages)} required stages did not produce a result — "
            + "; ".join(completeness.missing)
            + f"\n{_describe_completeness(completeness)}"
        )
    else:
        reporter.ok(
            f"done in {human_duration(total)} — all {len(completeness.stages)} stages ran — "
            f"{', '.join(produced)} in {out_dir}"
        )

    write_run_log(
        out_dir / "run.log",
        reporter_lines="",
        extra={
            **{k: str(v) for k, v in run_meta.items()},
            "total wall clock": human_duration(total),
            "reference entries": "; ".join(
                f"{k}={v.n_sequences}@{v.mean_length_aa:.0f}aa"
                for k, v in sorted(database.entries.items())
            ),
        },
    )
    reporter.write_log(out_dir / "run.log.tail")
    _merge_logs(out_dir)

    # ---- stdout ----------------------------------------------------------- #
    style = Style(sys.stdout.isatty() and not os.environ.get("NO_COLOR"))
    if not print_summary:
        pass
    elif args.json:
        print(json.dumps(results_json, indent=2))
    elif args.compact:
        print(
            render_compact(
                sample=sample_input.sample,
                tally=result,
                selected=selected_results,
                style=style,
            )
        )
    else:
        print(
            render_summary(
                sample=sample_input.sample,
                tally=result,
                selected=selected_results,
                qc=qc,
                profile=profile,
                diagnostics=diagnostics,
                run_meta=run_meta,
                style=style,
                similarity=similarity,
            )
        )

    if pdf_error is not None:
        return 1
    if not completeness.complete:
        return 2
    return 1 if any(d.level == "error" for d in diagnostics) else 0


#: Flags whose values are absolute paths that change between runs; showing them
#: in the report adds noise and breaks reproducibility comparisons.
_PATH_FLAGS: Final = frozenset({"--db", "--query", "--out"})


def _flags_for_display(command: Sequence[str]) -> str:
    """Render the DIAMOND command with per-run paths elided."""
    if not command:
        return "n/a"
    parts: list[str] = []
    skip_next = False
    for token in command[1:]:  # drop the absolute executable path
        if skip_next:
            skip_next = False
            continue
        if token in _PATH_FLAGS:
            skip_next = True
            continue
        parts.append(token)
    return " ".join(parts)


def _merge_logs(out_dir: Path) -> None:
    """Append the captured progress lines to run.log and drop the temp file."""
    tail = out_dir / "run.log.tail"
    if not tail.is_file():
        return
    log = out_dir / "run.log"
    with log.open("a", encoding="utf-8") as fh:
        fh.write(tail.read_text(encoding="utf-8"))
    tail.unlink(missing_ok=True)


# --------------------------------------------------------------------------- #
# entry point
# --------------------------------------------------------------------------- #


def cmd_fmt(args: argparse.Namespace) -> int:
    """`openbiota fmt` — FMT donor and donor-set matching (BUILD_SPEC_v06.0)."""
    from openbiota.fmt.cli import run as _run_fmt

    return _run_fmt(args)


def cmd_biofilm(args: argparse.Namespace) -> int:
    """`openbiota biofilm` — the biofilm module (BUILD_SPEC_v08.0 §15)."""
    from openbiota.biofilm.cli import run as _run_biofilm

    return _run_biofilm(args)


def cmd_mycobiome(args: argparse.Namespace) -> int:
    """`openbiota mycobiome` — the fungal module (BUILD_SPEC_v08.2)."""
    from openbiota.mycobiome.cli import run as _run_mycobiome

    return _run_mycobiome(args)


def main(argv: Sequence[str] | None = None) -> int:
    raw = list(sys.argv[1:] if argv is None else argv)
    # `openbiota` with no subcommand, or with only run-flags, means `openbiota run`.
    if not raw or (raw[0] not in SUBCOMMANDS and raw[0] not in {"-h", "--help", "--version"}):
        raw = ["run", *raw]

    parser = build_parser()
    args = parser.parse_args(raw)
    command = args.command or "run"

    handlers = {
        "run": cmd_run,
        "panels": cmd_panels,
        "build-db": cmd_build_db,
        "doctor": cmd_doctor,
        "cohort": cmd_cohort,
        "validate": cmd_validate,
        "depth-check": cmd_depth_check,
        "taxonomic-cohort": cmd_taxonomic_cohort,
        "profiles": cmd_profiles,
        "build-host-index": cmd_build_host_index,
        "run-all": cmd_run_all,
        "prune": cmd_prune,
        "validate-profiles": cmd_validate_profiles,
        "age-model": cmd_age_model,
        "fmt": cmd_fmt,
        "biofilm": cmd_biofilm,
        "mycobiome": cmd_mycobiome,
    }
    try:
        return handlers[command](args)
    except OpenBiotaError as exc:
        print(f"\nopenbiota: {exc}", file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        print("\nopenbiota: interrupted. Cached stages will be reused on the next run.", file=sys.stderr)
        return 130
    except subprocess.CalledProcessError as exc:  # pragma: no cover — defensive
        print(f"\nopenbiota: external command failed: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
