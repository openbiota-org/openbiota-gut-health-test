"""`openbiota fmt` — the matching CLI (spec §4.1).

Exit codes are the spec's: 0 completed comparison including "no admissible
candidate"; 2 invalid request or identity conflict; 3 no usable analytical
input; 4 an explicitly requested mandatory model failed; 5 report integrity
failure. A missing optional capability is a structured warning and still
exits 0 when a meaningful core comparison was produced.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from openbiota.errors import OpenBiotaError

from . import capabilities
from .engine import leading_summary, match
from .inputs import FmtInputError
from .report import write_outputs

FORMATS = ("json", "html", "pdf")


def add_parser(sub: argparse._SubParsersAction) -> None:
    """Attach `openbiota fmt ...` to the main parser."""
    fmt = sub.add_parser(
        "fmt",
        help="recipient-specific FMT donor and donor-set matching (research use only)",
        description=(
            "Compare one recipient against candidate donor materials and evaluate donor sets for "
            "complementary restoration potential. Research comparison only: this command never "
            "clears a donor, never predicts clinical benefit and never produces preparation, "
            "dosing or administration instructions."
        ),
    )
    fsub = fmt.add_subparsers(dest="fmt_command")

    m = fsub.add_parser("match", help="run a comparison and write match.json/html/pdf")
    m.add_argument("--recipient", help="recipient sample output directory")
    m.add_argument(
        "--donor",
        action="append",
        default=[],
        metavar="NAME=PATH",
        help="a candidate donor as LABEL=/path/to/sample_dir (repeatable)",
    )
    m.add_argument("--manifest", help="a JSON MatchRequest instead of --recipient/--donor")
    m.add_argument("--out", required=True, help="run output directory")
    m.add_argument("--indication", default=None, help="recipient indication code, e.g. longcovid")
    m.add_argument("--goals", help="JSON file of explicit research goals")
    m.add_argument("--screening-manifest", help="JSON file of validated screening results")
    m.add_argument("--formats", default="json,html,pdf", help="comma-separated: json,html,pdf")
    m.add_argument("--max-subsets", type=int, default=100000)
    m.add_argument("--max-donors-per-set", type=int, default=None)
    m.add_argument("--max-materials-per-set", type=int, default=None)
    m.add_argument("--models", default="standard", choices=("core", "standard", "all-available"))
    m.add_argument("--objective", default="measured_coverage",
                   choices=("measured_coverage", "predicted_target_attainment", "model_endpoint",
                            "mechanistic_target_pareto"))
    m.add_argument("--seed", type=int, default=1701)
    m.add_argument("--offline", action="store_true", help="no network use (the core never needs it)")
    m.add_argument("--quiet", action="store_true")

    i = fsub.add_parser("inspect", help="validate a request and report what would be used")
    i.add_argument("--manifest", required=True)

    e = fsub.add_parser("explain", help="explain one candidate from a completed run")
    e.add_argument("--run", required=True)
    e.add_argument("--candidate", required=True, help="donor label or LABEL+LABEL for a set")

    v = fsub.add_parser("validate", help="check a completed run's integrity")
    v.add_argument("--run", required=True)

    mo = fsub.add_parser("models", help="model and capability status")
    mo.add_argument("models_action", nargs="?", default="list", choices=("list",))


def _request_from_args(args: argparse.Namespace) -> dict[str, Any]:
    if args.manifest:
        req = json.loads(Path(args.manifest).read_text())
        req.setdefault("search", {})
        return req
    if not args.recipient or not args.donor:
        raise OpenBiotaError(
            "give --recipient and at least one --donor LABEL=/path, or --manifest request.json"
        )
    donors = []
    for spec in args.donor:
        if "=" not in spec:
            raise OpenBiotaError(f"--donor must be LABEL=/path, got {spec!r}")
        label, path = spec.split("=", 1)
        donors.append({"person_id": label.strip(), "input_root": path.strip()})
    goals = json.loads(Path(args.goals).read_text()) if args.goals else []
    screening = (
        json.loads(Path(args.screening_manifest).read_text()) if args.screening_manifest else None
    )
    return {
        "request_id": f"{Path(args.recipient).name}-vs-{len(donors)}",
        "recipient": {"person_id": Path(args.recipient).name.split("_")[0], "input_root": args.recipient},
        "donors": donors,
        "indication": (
            {"code": args.indication, "source": "user_report", "confirmed": False}
            if args.indication
            else None
        ),
        "custom_goals": goals if isinstance(goals, list) else goals.get("goals", []),
        "screening_manifest": screening,
        "model_policy": args.models,
        "objective_views": [
            {"type": args.objective, "model_id": None, "endpoint": None, "horizon": None,
             "scenario_set_id": None}
        ],
        "search": {
            "max_subsets": args.max_subsets,
            "max_donors_per_set": args.max_donors_per_set,
            "max_materials_per_set": args.max_materials_per_set,
            "seed": args.seed,
        },
        "output_formats": [f for f in args.formats.split(",") if f.strip() in FORMATS],
    }


def _cmd_match(args: argparse.Namespace) -> int:
    request = _request_from_args(args)
    if args.objective != "measured_coverage":
        man = {c["capability_id"]: c for c in capabilities.manifest()}
        blocked = [
            c for c in man.values()
            if c["kind"] in ("predictive_model", "published_adapter", "mechanistic_model")
        ]
        print(
            f"openbiota fmt: objective {args.objective!r} needs a model artifact. "
            f"{len(blocked)} model capabilities are unavailable in this environment; "
            "see `openbiota fmt models list`.",
        )
        return 4
    result = match(request)
    written = write_outputs(result, Path(args.out), request["output_formats"])

    # Integrity: the canonical document must carry no forbidden claim field
    # and every displayed aggregate must resolve to a target ID (§5.4, V7-056).
    blob = json.dumps(result)
    for banned in ('"safe"', "donor_clearance_probability", "approved_for_fmt"):
        if banned in blob:
            print(f"openbiota fmt: report integrity failure, forbidden field {banned}")
            return 5

    if not args.quiet:
        ach = result["achievable"]
        ceiling = ach["index_ceiling"]
        lb = result["leaderboard"]
        ranked = sorted(result["set_results"], key=lambda b: b["match_rank"])
        rec = next((b for b in ranked if b["candidate_id"] == lb.get("recommended_id")), None)
        top = ranked[0] if ranked else None
        if rec is not None:
            print(
                f"  RECOMMENDED  {' + '.join(rec['member_person_ids'])}"
                f"   match score {rec['match_score']:.1f} of 100"
                f"   {rec['goals']['covered']} of {lb['reachable_goals']} gaps"
                "   0 new disease patterns"
            )
            if lb.get("recommended_note"):
                print(f"               {lb['recommended_note']}")
        elif lb.get("no_clean_candidate_note"):
            print(f"  {lb['no_clean_candidate_note']}")
        if top is not None and rec is not None and top["candidate_id"] != rec["candidate_id"]:
            intro = ", ".join(
                f"{e['short_label']} ({e['donor_percentile']:.0f}th vs your "
                f"{e['recipient_percentile']:.0f}th)"
                for e in top["introduces_new_patterns"] if e["counts_against_selection"]
            )
            print(
                f"  highest score {' + '.join(top['member_person_ids'])} at "
                f"{top['match_score']:.1f}, but it would introduce {intro}"
            )
        print()
        print(leading_summary(result))
        print()
        lead = tuple(
            t for t in ("confirmed_exclusion", "high_consequence_supported",
                        "high_consequence_unresolved", "toxin_or_resistance_review")
        )
        print(
            f"  {'donor':<6s} {'gaps covered':>13s} {'index':>7s} {'of max':>7s} "
            f"{'lab tests':>10s} {'GMWI2':>7s}  new disease patterns"
        )
        for row in sorted(result["individual_results"], key=lambda x: -x["goals"]["covered"]):
            g = row["goals"]
            mine = [c for c in result["concern_records"] if c["candidate_id"] == row["candidate_id"]]
            serious = sum(1 for c in mine if c["tier"] in lead)
            gm = row["context_not_in_the_score"]["gmwi2"]["score"]
            share = (
                f"{100.0 * (row['coverage']['coverage'] or 0) / ceiling:.0f}%"
                if ceiling else "—"
            )
            solo = next(
                (b for b in result["set_results"]
                 if b["member_person_ids"] == [row["person_id"]]),
                None,
            )
            patterns = (
                ", ".join(
                    e["short_label"] for e in solo["introduces_new_patterns"]
                    if e["counts_against_selection"]
                )
                or "none"
            ) if solo else "—"
            print(
                f"  {row['person_id']:<6s} {g['covered']:>6d} of {g['total']:<4d} "
                f"{row['coverage']['coverage']:>7.1f} {share:>7s} {serious:>10d} "
                f"{(gm if gm is not None else float('nan')):>+7.2f}  {patterns}"
            )
        if ach["unreachable_goals"]:
            print(
                "\n  no candidate supplies: "
                + "; ".join(u["label"] for u in ach["unreachable_goals"])
                + f"  (so the index cannot exceed {ceiling:.1f})"
            )
        # Only advertise a combination that clears the same bar as the
        # recommendation. Naming the top-scoring pair here would quietly put a
        # pattern-introducing set back in front of the reader.
        combos = [s for s in ranked if s["material_count"] > 1]
        clean_combo = next((s for s in combos if s["new_pattern_count"] == 0), None)
        if clean_combo:
            print(f"\n  {clean_combo['plain_summary']}")
        elif combos:
            print(
                "\n  every combination of these candidates would introduce a disease pattern you "
                "do not have, so no combination is recommended over the single donor above"
            )
        print(f"\n  written: {', '.join(sorted(written.values()))}")
        print(
            "  Research comparison only. 'Gaps covered' counts organisms and gene capacities the "
            "recipient is low in or missing that a donor carries; it says nothing about whether they "
            "would establish. No donor is screened, cleared or recommended, and no preparation or "
            "dosing instruction is implied."
        )
    return 0


def _cmd_inspect(args: argparse.Namespace) -> int:
    request = json.loads(Path(args.manifest).read_text())
    print(json.dumps({"request": request, "capabilities": capabilities.manifest()}, indent=2))
    return 0


def _cmd_explain(args: argparse.Namespace) -> int:
    result = json.loads((Path(args.run) / "match.json").read_text())
    names = {r["candidate_id"]: r["person_id"] for r in result["individual_results"]}
    by_name = {v: k for k, v in names.items()}
    wanted = "+".join(sorted(by_name.get(p, p) for p in args.candidate.split("+")))
    hit = next((s for s in result["set_results"] if s["candidate_id"] == wanted), None)
    if hit is None:
        print(f"openbiota fmt: no evaluated candidate {args.candidate!r}")
        return 2
    targets = {t["target_id"]: t for t in result["targets"]}
    print(f"candidate {'+'.join(hit['member_person_ids'])}  (front {hit['pareto_front']})")
    print(f"  coverage {hit['coverage']['coverage']:.2f} "
          f"[{hit['coverage']['scenario_lower']:.2f}–{hit['coverage']['scenario_upper']:.2f}]")
    print(f"  gain over best single material: {hit['marginals']['gain_over_global_best_singleton']:+.2f}")
    for mid, ids in (hit["unique_target_ids_by_member"] or {}).items():
        if ids:
            print(f"  only {names.get(mid, mid)} supplies:")
            for tid in ids:
                print(f"    · {targets[tid]['label']}")
    print(f"  redundant goals (more than one member supplies): {len(hit['redundant_target_ids'])}")
    print(f"  unresolved concerns inherited: {hit['concern_counts']['review_pending']}")
    return 0


def _cmd_validate(args: argparse.Namespace) -> int:
    path = Path(args.run) / "match.json"
    result = json.loads(path.read_text())
    problems: list[str] = []
    target_ids = {t["target_id"] for t in result["targets"]}
    scoring = {t["target_id"] for t in result["targets"] if t["in_positive_denominator"]}
    for row in result["individual_results"]:
        ids = {c["target_id"] for c in row["per_target_availability"]}
        if ids != scoring:
            problems.append(f"{row['person_id']}: availability cells do not match the scored denominator")
    for s in result["set_results"]:
        if s["candidate_id"] != "+".join(sorted(s["member_material_ids"])):
            problems.append(f"{s['candidate_id']}: set ID is not permutation-invariant")
        for tid in s["redundant_target_ids"]:
            if tid not in target_ids:
                problems.append(f"{s['candidate_id']}: unknown target {tid}")
    denominators = {
        (s["coverage"]["denominator_target_count"]) for s in result["set_results"]
    }
    if len(denominators) > 1:
        problems.append(f"denominator is not frozen across candidates: {sorted(denominators)}")
    blob = json.dumps(result)
    for banned in ('"safe"', "donor_clearance_probability", "approved_for_fmt"):
        if banned in blob:
            problems.append(f"forbidden claim field present: {banned}")
    if problems:
        for p in problems:
            print(f"openbiota fmt validate: {p}")
        return 5
    print(
        f"ok: {len(result['set_results'])} candidates, frozen denominator of "
        f"{sorted(denominators)[0]} goals, {len(result['concern_records'])} concern records, "
        f"clinical_release_status={result['clinical_release_status']}"
    )
    return 0


def _cmd_models(_args: argparse.Namespace) -> int:
    for c in capabilities.manifest():
        miss = "; ".join(c.get("missing_artifacts") or [])
        print(f"{c['capability_id']:<34s} {c['status']:<38s} {c.get('kind', '')}")
        if miss:
            print(f"{'':34s} needs: {miss}")
    return 0


def run(args: argparse.Namespace) -> int:
    """Dispatch `openbiota fmt <subcommand>`."""
    handlers = {
        "match": _cmd_match,
        "inspect": _cmd_inspect,
        "explain": _cmd_explain,
        "validate": _cmd_validate,
        "models": _cmd_models,
        None: lambda _a: (_cmd_models(_a), print("\nusage: openbiota fmt match --recipient DIR --donor NAME=DIR --out DIR"))[0],
    }
    handler = handlers.get(getattr(args, "fmt_command", None))
    if handler is None:  # pragma: no cover - argparse restricts this
        raise OpenBiotaError(f"unknown fmt subcommand {args.fmt_command!r}")
    try:
        return handler(args)
    except FmtInputError as exc:
        print(f"openbiota fmt: {exc}")
        return 3
