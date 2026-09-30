"""Did every stage of the report actually run?

A run used to be able to finish "successfully" with a stage quietly
skipped: every stage is wrapped so that a missing tool or a failure warns
and moves on, which is right for keeping the rest of the report alive but
wrong as the last word. One sample's report said "strain typing did not
run" on page one; another's simulation was never solved because the run
was in a cached mode. Both runs exited 0.

This module is the last word. It reads the finished `results.json` - the
same file the report is rendered from - and asks each required stage the
one question that matters: did you produce your result? The answer goes
into `run.completeness` in the file, the missing stages are printed in red
at the end of the run, and the run exits non-zero if anything is missing.

A stage is checked by its *output*, never by whether its code was reached:
a stage that ran and produced nothing is a stage that did not run.

To add a stage: append a `Stage` to `REQUIRED`. The check receives the
whole results object and returns ``(ran, detail)``; keep the detail short,
it is printed on one line.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from typing import Any, Final

Check = Callable[[Mapping[str, Any]], tuple[bool, str]]


@dataclass(frozen=True)
class Stage:
    key: str
    label: str
    check: Check


def _d(results: Mapping[str, Any], *path: str) -> Any:
    """Walk a path of keys, returning None the moment one is missing."""
    node: Any = results
    for key in path:
        if not isinstance(node, Mapping):
            return None
        node = node.get(key)
    return node


def _sequencing(r: Mapping[str, Any]) -> tuple[bool, str]:
    gates = _d(r, "sequencing_quality", "gates", "gates") or []
    return bool(gates), f"{len(gates)} quality gates"


def _functional(r: Mapping[str, Any]) -> tuple[bool, str]:
    panels = r.get("panels") or []
    compared = _d(r, "reference_comparison", "panels") or {}
    ok = bool(panels) and bool(compared)
    return ok, f"{len(panels)} panels, {len(compared)} placed against the reference"


def _taxonomy(r: Mapping[str, Any]) -> tuple[bool, str]:
    ok = bool(_d(r, "community_profile", "available"))
    return ok, str(_d(r, "community_profile", "basis") or "no community profile")


def _extended(r: Mapping[str, Any]) -> tuple[bool, str]:
    ext = r.get("extended_catalogue") or {}
    ok = ext.get("engine") == "metaphlan4" and bool(ext.get("kingdom_counts") or ext.get("n_genera_detected"))
    return ok, f"{ext.get('database', 'no MetaPhlAn 4 catalogue')}"


def _genome(r: Mapping[str, Any]) -> tuple[bool, str]:
    g = r.get("genome_profile") or {}
    return g.get("status") == "resolved", f"{g.get('n_hits', 0)} species clusters"


def _strain(r: Mapping[str, Any]) -> tuple[bool, str]:
    s = r.get("strain_resolution") or {}
    return s.get("status") == "resolved", (
        f"{s.get('organisms_resolved', 0)} organisms resolved"
        if s.get("status") == "resolved" else f"status {s.get('status', 'absent')}"
    )


def _pathogens(r: Mapping[str, Any]) -> tuple[bool, str]:
    p = r.get("pathogens") or {}
    return p.get("analysis_status") == "completed", f"status {p.get('analysis_status', 'absent')}"


def _age(r: Mapping[str, Any]) -> tuple[bool, str]:
    a = r.get("microbiome_age") or {}
    # Abstaining on an out-of-distribution sample is the model doing its
    # job; only "not computable" and absence are failures to run.
    return a.get("status") in ("scored", "abstained_ood"), f"status {a.get('status', 'absent')}"


def _similarity(r: Mapping[str, Any]) -> tuple[bool, str]:
    ranked = _d(r, "profile_similarity", "ranked") or []
    return bool(ranked), f"{len(ranked)} published patterns scored"


def _mycobiome(r: Mapping[str, Any]) -> tuple[bool, str]:
    m = r.get("mycobiome") or {}
    return m.get("analysis_status") == "complete", f"status {m.get('analysis_status', 'absent')}"


def _biofilm(r: Mapping[str, Any]) -> tuple[bool, str]:
    b = r.get("biofilm") or {}
    scored = [c for c in (b.get("cards") or []) if c.get("reference_percentile") is not None]
    return b.get("profile_status") == "available" and bool(scored), f"{len(scored)} cards scored"


def _findings(r: Mapping[str, Any]) -> tuple[bool, str]:
    f = r.get("findings_and_evidence") or {}
    ok = isinstance(f.get("triggers"), list) and bool(f.get("triggers")) and "findings" in f
    return ok, f"{len(f.get('triggers') or [])} triggers, {len(f.get('findings') or [])} findings"


def _inventory(r: Mapping[str, Any]) -> tuple[bool, str]:
    inv = r.get("organism_inventory") or {}
    return bool(inv.get("n_organisms")), f"{inv.get('n_organisms', 0)} organisms"


def _extension(r: Mapping[str, Any]) -> tuple[bool, str]:
    scores = _d(r, "extension", "views", "reading_scores") or {}
    return bool(scores), f"{len(scores)} readings scored"


def _simulation(r: Mapping[str, Any]) -> tuple[bool, str]:
    state = _d(r, "extension", "views", "simulation", "scenarios", "execution_state")
    return state == "solved_model", f"execution state {state or 'absent'}"


def _assembly(r: Mapping[str, Any]) -> tuple[bool, str]:
    a = r.get("assembly") or {}
    st = a.get("status")
    # Not triggered is a complete answer: nothing needed assembling.
    ok = st in ("completed", "not_triggered")
    if st == "completed":
        return ok, f"{a['assembly']['n_contigs']:,} contigs, N50 {a['assembly']['n50']:,}"
    return ok, f"status {st or 'absent'}: {a.get('reason', '')}".strip(": ")


def _bacdive(r: Mapping[str, Any]) -> tuple[bool, str]:
    b = r.get("bacdive") or {}
    ok = "species" in b and b.get("status") != "failed"
    return ok, (f"{b.get('n_found', 0)} of {b.get('n_species_queried', 0)} named species found" if ok
                else f"status {b.get('status', 'absent')}")


def _strain_analysis(r: Mapping[str, Any]) -> tuple[bool, str]:
    s = r.get("strain_analysis") or {}
    ok = s.get("status") == "completed"
    return ok, (f"{s.get('n_compared', 0)} placed, {s.get('n_unresolved', 0)} unresolved" if ok
                else f"status {s.get('status', 'absent')}: {s.get('reason', '')}".strip(": "))


def _confirmation(r: Mapping[str, Any]) -> tuple[bool, str]:
    c = _d(r, "detection", "confirmation") or {}
    ok = c.get("status") == "completed"
    vs = c.get("verdicts") or []
    return ok, (f"{c.get('n_tested', 0)} tested, {sum(1 for v in vs if v['status'] == 'pending')} pending a genome"
                if ok else f"status {c.get('status', 'absent')}")


def _lane(lane_id: str) -> Check:
    def check(r: Mapping[str, Any]) -> tuple[bool, str]:
        st = _d(r, "detection", "lane_status", lane_id) or {}
        n = st.get("n_species")
        return st.get("state") == "ran", (
            f"{n} species-level observations" if st.get("state") == "ran"
            else f"{st.get('state', 'absent')}: {st.get('reason', '')}".strip(": ")
        )
    return check


REQUIRED: Final[tuple[Stage, ...]] = (
    Stage("sequencing_quality", "sequencing quality", _sequencing),
    Stage("functional_search", "functional search and reference placement", _functional),
    Stage("taxonomy", "taxonomic profile", _taxonomy),
    Stage("extended_catalogue", "extended catalogue (MetaPhlAn 4)", _extended),
    Stage("genome_profile", "whole-genome profile", _genome),
    Stage("strain_typing", "strain typing (marker consensus)", _strain),
    Stage("pathogens", "pathogen screen", _pathogens),
    Stage("age", "estimated age of biota", _age),
    Stage("profile_similarity", "resemblance to published patterns", _similarity),
    Stage("mycobiome", "mycobiome", _mycobiome),
    Stage("biofilm", "biofilm-related potential", _biofilm),
    Stage("findings", "organism findings and evidence", _findings),
    Stage("organism_inventory", "organism inventory", _inventory),
    Stage("reading_scores", "reading scores", _extension),
    Stage("simulation", "model-assisted simulation", _simulation),
    # Spec 0.8.4: the expanded detection lanes. Each judged by whether it
    # produced observations; a missing database is MISSING, never silent.
    Stage("lane_jan26", "marker lane (MetaPhlAn 4.2.6 / Jan26)", _lane("metaphlan_jan26")),
    Stage("lane_globdb", "broad discovery (sylph / GlobDB r232)", _lane("sylph_globdb")),
    Stage("lane_motus", "universal-marker lane (mOTUs 4.1)", _lane("motus4")),
    Stage("lane_kraken", "gut rescue classifier (Kraken2 + Bracken / UHGG)", _lane("kraken_uhgg")),
    Stage("lane_singlem", "unrepresented-lineage rescue (SingleM / GlobDB)", _lane("singlem_globdb")),
    Stage("confirmation", "competitive confirmation of new and single-method calls", _confirmation),
    Stage("strain_analysis", "comparative strain analysis (StrainPhlAn 4 / Jan26)", _strain_analysis),
    Stage("bacdive", "BacDive enrichment", _bacdive),
    Stage("assembly", "targeted assembly of unresolved lineages", _assembly),
)


@dataclass
class Completeness:
    stages: list[dict[str, Any]]

    @property
    def complete(self) -> bool:
        return all(s["ran"] for s in self.stages)

    @property
    def missing(self) -> list[str]:
        return [s["label"] for s in self.stages if not s["ran"]]

    def to_json(self) -> dict[str, Any]:
        return {
            "complete": self.complete,
            "n_required": len(self.stages),
            "n_ran": sum(1 for s in self.stages if s["ran"]),
            "missing": self.missing,
            "stages": list(self.stages),
            "meaning": (
                "Every required stage of the report, and whether it produced its result in "
                "this run. A stage is judged by its output in this file, not by whether it "
                "was attempted. A run with anything missing exits non-zero."
            ),
        }


def audit(results: Mapping[str, Any]) -> Completeness:
    """Ask every required stage whether it produced its result."""
    rows: list[dict[str, Any]] = []
    for stage in REQUIRED:
        try:
            ran, detail = stage.check(results)
        except Exception as exc:  # noqa: BLE001 - a broken check is a failed stage
            ran, detail = False, f"check failed: {type(exc).__name__}: {exc}"
        rows.append({"key": stage.key, "label": stage.label, "ran": bool(ran), "detail": detail})
    return Completeness(rows)


def describe(report: Completeness) -> str:
    """One line per stage, for the end of the run."""
    width = max(len(s["label"]) for s in report.stages)
    return "\n".join(
        f"  {'ran    ' if s['ran'] else 'MISSING'}  {s['label']:<{width}}  {s['detail']}"
        for s in report.stages
    )
