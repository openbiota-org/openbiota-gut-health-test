"""`openbiota biofilm` — the biofilm module's command surface.

Each verb has the semantics section 15 of the spec assigns it. Two are
worth stating because they are easy to get subtly wrong:

``compare`` prints the recipient and each donor side by side and stops
there. It does not compute a pooled result, a combined inventory
percentage or a safety verdict, because relative abundances are
compositional and pooled material cannot be inferred by adding them.

``datasets audit`` reports that zero of the fourteen registered datasets
are training-eligible. That is the correct answer today and printing it
plainly is the point: a future change that starts training on biopsy 16S
should have to change this number visibly.

Exit codes follow the spec: 0 for a completed analysis even when it has
explicit partial results, nonzero only for corrupt input or execution
failure. A missing optional assay is not a crash.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from . import datasets as ds_mod
from . import engine, interventions, reference, registry, scoring
from . import sources as src


class BiofilmCliError(RuntimeError):
    """A usage or input error worth a nonzero exit."""


def _resolve_sample(args: argparse.Namespace, sample: str) -> Path:
    """Accept either a sample ID or a path to a results directory."""
    candidate = Path(sample)
    if candidate.is_dir() and (candidate / "results.json").exists():
        return candidate
    guess = Path(args.results_dir) / sample
    if (guess / "results.json").exists():
        return guess
    # Tolerate an ID prefix, so `--sample K1` finds `SAMPLE4_A04`.
    matches = sorted(
        p
        for p in Path(args.results_dir).glob(f"{sample}*")
        if (p / "results.json").exists()
    )
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise BiofilmCliError(
            f"{sample!r} matches {len(matches)} results directories: "
            + ", ".join(p.name for p in matches)
        )
    raise BiofilmCliError(f"no results found for {sample!r} under {args.results_dir}")


def _load(path: Path) -> dict[str, Any]:
    try:
        return json.loads((path / "results.json").read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise BiofilmCliError(f"cannot read {path / 'results.json'}: {exc}") from exc


def _analyze(args: argparse.Namespace, sample: str) -> dict[str, Any]:
    path = _resolve_sample(args, sample)
    results = _load(path)
    cached = results.get("biofilm")
    if isinstance(cached, dict) and cached.get("cards"):
        return cached
    return engine.analyze(path, results, cohort_path=str(args.cohort))


def _emit(payload: dict[str, Any], args: argparse.Namespace, text: list[str]) -> int:
    if getattr(args, "output_format", "text") == "json":
        print(json.dumps(payload, indent=2))
    else:
        print("\n".join(text))
    return 0


def _card_lines(doc: dict[str, Any]) -> list[str]:
    out = []
    for card in doc.get("cards", []):
        pct = card.get("reference_percentile")
        value = f"{pct:5.1f}" if pct is not None else "    -"
        out.append(
            f"  {card['card']:4s} {card['label'][:44]:44s} {value}  "
            f"{card['band'] if pct is not None else card['status']}"
        )
    return out


def cmd_inventory(args: argparse.Namespace) -> int:
    """Enumerate inputs and reference compatibility before expensive work."""
    path = _resolve_sample(args, args.sample)
    mpa3 = path / "taxonomy" / reference.MPA3_PROFILE
    mpa4 = path / "taxonomy" / reference.MPA4_PROFILE
    try:
        cohort = reference.load_reference(str(args.cohort))
        cohort_state: dict[str, Any] = cohort.provenance()
    except reference.ReferenceUnavailable as exc:
        cohort_state = {"available": False, "reason": str(exc)}

    payload = {
        "sample": path.name,
        "results_json": (path / "results.json").exists(),
        "profiles": {
            "metaphlan3_reference_compatible": {
                "path": str(mpa3),
                "present": mpa3.exists(),
                "used_for": "the H-C and P-E proxies",
            },
            "metaphlan4_sgb": {
                "path": str(mpa4),
                "present": mpa4.exists(),
                "used_for": (
                    "organism inventory elsewhere in the report; NOT ranked "
                    "against the MetaPhlAn 3 reference, because they are "
                    "different measurements"
                ),
            },
        },
        "reference": cohort_state,
        "registry": {
            "modules": len(registry.MODULES),
            "sources": len(src.SOURCES),
            "datasets": len(ds_mod.DATASETS),
            "interventions": len(interventions.INTERVENTIONS),
        },
    }
    text = [
        f"sample                 {path.name}",
        f"MetaPhlAn 3 profile    {'present' if mpa3.exists() else 'MISSING'}"
        "   (the reference-compatible lane)",
        f"MetaPhlAn 4 profile    {'present' if mpa4.exists() else 'missing'}"
        "   (not used for these proxies)",
        f"reference cohort       {cohort_state.get('n_participants', 'unavailable')}"
        " independent adults",
        f"registry               {len(registry.MODULES)} modules, "
        f"{len(src.SOURCES)} sources, {len(ds_mod.DATASETS)} datasets, "
        f"{len(interventions.INTERVENTIONS)} interventions",
    ]
    return _emit(payload, args, text)


def cmd_analyze(args: argparse.Namespace) -> int:
    path = _resolve_sample(args, args.sample)
    doc = engine.analyze(path, _load(path), cohort_path=str(args.cohort))
    text = [f"{doc['sample_id']} — {doc['measurement_label']}", ""]
    text += _card_lines(doc)
    text.append("")
    text.append(f"  {doc['headline_note']}")
    if doc.get("out_of_domain"):
        text.append(f"  out of domain: {doc['out_of_domain']}")
    return _emit(doc, args, text)


def cmd_explain(args: argparse.Namespace) -> int:
    """Show the evidence behind a result without rerunning alignment."""
    doc = _analyze(args, args.sample)
    text = [f"{doc['sample_id']} — how each number was produced", ""]
    for card in doc.get("cards", []):
        text.append(f"{card['card']} · {card['label']}")
        pct = card.get("reference_percentile")
        if pct is None:
            text.append(f"    not scored: {', '.join(card.get('reason_codes') or [])}")
            for omission in card.get("omissions", []):
                text.append(f"      - {omission}")
            text.append("")
            continue
        text.append(
            f"    percentile {pct:.2f} against {card['reference_n']:,} adults "
            f"(calibration {card['calibration_id']})"
        )
        for name, value in (card.get("feature_values") or {}).items():
            fp = (card.get("feature_percentiles") or {}).get(name)
            d = (card.get("feature_detail") or {}).get(name, {})
            mark = "" if d.get("detected") else "   [not detected above assay limit]"
            text.append(f"      {name:18s} {value:9.4f}%  rank {fp:6.2f}{mark}")
        driver = card.get("largest_driver") or {}
        if driver:
            text.append(f"    largest driver: {driver.get('note', '')}")
        text.append(f"    not: {card.get('what_it_is_not', '')}")
        text.append("")
    if getattr(args, "include_candidate_evidence", False):
        text.append("Registered mechanisms not currently scoring:")
        for m in registry.MODULES.values():
            if not m.eligible_for_axis:
                text.append(f"    {m.id}  {m.label[:48]:48s} {m.ineligibility_reason}")
    return _emit(doc, args, text)


def cmd_rank_actions(args: argparse.Namespace) -> int:
    """Rank matched evidence from an existing result; no sequencing rerun."""
    doc = _analyze(args, args.sample)
    actions = doc.get("intervention_candidates", [])
    text = [f"{doc['sample_id']} — {len(actions)} matched evidence records", ""]
    section_titles = {
        "support_protective_ecology": "Support protective ecology",
        "reduce_matched_concerning_mechanism": "Reduce a matched concerning mechanism",
        "clarify_or_treat_established_condition": (
            "Clarify or treat an independently established condition"
        ),
        "investigational_precision": "Investigational precision approaches",
    }
    for key, title in section_titles.items():
        rows = [a for a in actions if a["section"] == key]
        if not rows:
            continue
        text.append(f"{title} ({len(rows)})")
        for a in rows:
            text.append(
                f"    {a['id']}  {a['compound'][:38]:38s} "
                f"{a['endpoint_plain'][:30]:30s} {a['review_label']}"
            )
        text.append("")
    text.append(
        "These are records to read, not instructions. Laboratory exposures are "
        "not doses."
    )
    return _emit({"sample_id": doc["sample_id"], "actions": actions}, args, text)


def cmd_compare(args: argparse.Namespace) -> int:
    """Recipient and donors side by side, with nothing pooled."""
    names = [args.recipient, *args.donors]
    docs = {}
    for name in names:
        path = _resolve_sample(args, name)
        docs[path.name] = engine.analyze(
            path, _load(path), cohort_path=str(args.cohort)
        )

    text = ["Biofilm-related potential, per person", ""]
    text.append(f"  {'person':22s} {'H-M':>6s} {'H-C':>6s} {'P-M':>6s} {'P-E':>6s}")

    def cell(cards: dict[str, Any], key: str) -> str:
        pct = cards.get(key, {}).get("reference_percentile")
        return f"{pct:6.1f}" if pct is not None else "     -"

    for name, doc in docs.items():
        cards = {c["card"]: c for c in doc["cards"]}
        text.append(
            f"  {name[:22]:22s} {cell(cards, 'H-M')} {cell(cards, 'H-C')} "
            f"{cell(cards, 'P-M')} {cell(cards, 'P-E')}"
        )
    text += [
        "",
        "These candidates differ in detected biofilm-associated mechanisms. The",
        "measurements do not establish whether their communities will form",
        "protective or harmful biofilms in this recipient. No validated numerical",
        "probability of transferring a biofilm-related disease can be calculated",
        "from these results.",
        "",
        "No pooled score is shown. Relative abundances are compositional, so the",
        "mixture cannot be inferred by adding percentages, and a protective",
        "reading in one person cannot offset a concerning reading in another.",
    ]
    payload = {
        "recipient": args.recipient,
        "donors": list(args.donors),
        "per_person": {k: v["summary"] for k, v in docs.items()},
        "pooled_result": None,
        "pooled_refused_because": (
            "Relative abundances are compositional and engraftment, competition, "
            "dose, viability and host conditions are unmeasured. A union of genes "
            "would predict only a possible combined inventory."
        ),
    }
    return _emit(payload, args, text)


def cmd_validate(args: argparse.Namespace) -> int:
    """Run implementation checks and report what remains unvalidated."""
    checks: list[tuple[str, bool, str]] = []

    def check(name: str, fn: Any, detail: str = "") -> None:
        try:
            fn()
            checks.append((name, True, detail))
        except Exception as exc:  # noqa: BLE001 - reporting is the point
            checks.append((name, False, f"{type(exc).__name__}: {exc}"))

    check("rank arithmetic fixtures (§8.6, §8.11)", scoring.self_test)
    check("reference construction invariants", reference.self_test)
    check("module registry builds", lambda: registry.census())
    check("intervention registry builds", lambda: interventions.census())
    check("dataset audit builds", lambda: ds_mod.audit())
    check(
        "every module source resolves",
        lambda: [src.resolve_all(m.primary_source_ids) for m in registry.MODULES.values()],
    )
    check(
        "opposing evidence closes",
        lambda: [
            interventions.with_opposing((i,)) for i in interventions.INTERVENTIONS
        ],
    )

    failed = [c for c in checks if not c[1]]
    text = [f"Implementation checks: {len(checks) - len(failed)}/{len(checks)} passed", ""]
    for name, ok, detail in checks:
        text.append(f"  {'PASS' if ok else 'FAIL'}  {name}")
        if not ok:
            text.append(f"        {detail}")
    text += [
        "",
        "Scientific validation that remains unavailable, and is not claimed:",
        "  - No validated stool-DNA classifier of protective versus harmful biofilm.",
        "  - No clinical outcome validation for either axis.",
        "  - No age-matched pediatric reference; the cohort is adults only.",
        "  - No mechanism module has a validated sequence panel, so H-M and P-M",
        "    are structurally unavailable rather than merely empty.",
        f"  - {ds_mod.audit()['n_training_eligible']} of "
        f"{len(ds_mod.DATASETS)} registered datasets are training-eligible.",
    ]
    payload = {
        "registry": args.registry,
        "suite": args.suite,
        "checks": [{"name": n, "passed": ok, "detail": d} for n, ok, d in checks],
        "n_passed": len(checks) - len(failed),
        "n_total": len(checks),
        "unavailable_validations": [
            "stool_dna_good_bad_classifier",
            "clinical_outcome_validation",
            "pediatric_reference",
            "mechanism_sequence_panels",
        ],
    }
    print(json.dumps(payload, indent=2)) if args.output_format == "json" else print(
        "\n".join(text)
    )
    return 1 if failed else 0


def cmd_datasets(args: argparse.Namespace) -> int:
    audit = ds_mod.audit()
    text = [
        f"Registered public datasets: {audit['n_datasets']}",
        f"Training-eligible:          {audit['n_training_eligible']}",
        "",
    ]
    for row in audit["datasets"]:
        text.append(
            f"  {row['id']}  {row['assay']:22s} {row['specimen']:20s} "
            f"{row['access_status']:16s} {row['label_mapping_status']}"
        )
        if row["blocking_reason"]:
            text.append(f"           blocked: {row['blocking_reason']}")
    text += ["", audit["note"]]
    return _emit(audit, args, text)


_HANDLERS = {
    "inventory": cmd_inventory,
    "analyze": cmd_analyze,
    "explain": cmd_explain,
    "rank-actions": cmd_rank_actions,
    "compare": cmd_compare,
    "validate": cmd_validate,
    "datasets": cmd_datasets,
}


def run(args: argparse.Namespace) -> int:
    command = getattr(args, "biofilm_cmd", None)
    if not command:
        print(
            "usage: openbiota biofilm {inventory,analyze,explain,compare,"
            "validate,rank-actions,datasets} ...",
        )
        return 2
    handler = _HANDLERS.get(command)
    if handler is None:
        print(f"openbiota biofilm: unknown subcommand {command!r}")
        return 2
    try:
        return handler(args)
    except BiofilmCliError as exc:
        print(f"openbiota biofilm: {exc}")
        return 2
