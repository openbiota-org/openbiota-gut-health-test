"""The mycobiome evidence registry and the MHS-E1 rule manifest.

Data lives in ``taxa/classification/mycobiome_evidence.yaml``. This module
resolves a detected fungus to its registry entry (by accepted name or
synonym: *Candida glabrata* and *Nakaseomyces glabratus* are one entry) and
evaluates every rule against a sample's findings, producing the
:class:`~openbiota.mycobiome.score.Contribution` records that feed the
scorer. Rules never look at read fractions or abundance for a presence or
genotype trigger; activation is 1 when the identity predicate is supported,
0 otherwise, and an unresolved higher-resolution claim enters only as an
admissible alternative.
"""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Any, Final

import yaml

from openbiota.mycobiome.score import Contribution, policy_lock

REGISTRY: Final = Path("taxa/classification/mycobiome_evidence.yaml")

EVIDENCE_LABELS: Final[dict[str, str]] = {
    "human_trial": "human trial",
    "human_association": "human association",
    "human_isolate_plus_experiment": "human isolate + experiment",
    "animal_experiment": "animal experiment",
    "in_vitro": "laboratory",
    "genome_inference": "genome inference",
    "insufficient_evidence": "insufficient evidence",
}

CATEGORY_LABELS: Final[dict[str, str]] = {
    "studied_probiotic_species": "species with studied probiotic strains",
    "commensal_opportunist": "common resident, opportunistic potential",
    "opportunist": "opportunistic yeast",
    "food_associated": "food-associated",
    "environmental_or_food": "environmental or food mould",
    "beneficial_preclinical": "preclinical beneficial evidence",
    "pathogen_alert_route": "pathogen (alert policy)",
    "unknown": "no health evidence",
}


@dataclass(frozen=True, slots=True)
class Claim:
    claim_id: str
    label: str
    direction: str
    endpoint: str
    host: str
    site: str
    design: str
    source: str
    limit: str
    detail: str = ""

    @property
    def label_words(self) -> str:
        return EVIDENCE_LABELS.get(self.label, self.label.replace("_", " "))


@dataclass(frozen=True, slots=True)
class TaxonEntry:
    name: str
    taxid: int | None
    rank: str
    synonyms: tuple[str, ...]
    category: str
    origin_default: str
    role: str
    evidence: tuple[Claim, ...]
    must_not_infer: str
    complex: str | None = None

    @property
    def category_words(self) -> str:
        return CATEGORY_LABELS.get(self.category, self.category.replace("_", " "))

    @property
    def favorable(self) -> tuple[Claim, ...]:
        return tuple(c for c in self.evidence if c.direction == "favorable")

    @property
    def concerning(self) -> tuple[Claim, ...]:
        return tuple(c for c in self.evidence if c.direction == "concerning")


@dataclass(frozen=True, slots=True)
class Rule:
    rule_id: str
    direction: str
    biological_claim_id: str
    source_ids: tuple[str, ...]
    required_identity: str
    taxon: str | None
    trigger_type: str
    cap: float
    bound: bool = True
    reference_accessions: tuple[str, ...] = ()
    reference_genotypes: tuple[str, ...] = ()
    co_present_bacteria: tuple[str, ...] = ()
    required_context: tuple[str, ...] = ()
    extrapolations: tuple[str, ...] = ()
    note: str = ""

    def to_json(self) -> dict[str, Any]:
        return {
            "rule_id": self.rule_id, "direction": self.direction, "biological_claim_id": self.biological_claim_id,
            "source_ids": list(self.source_ids), "required_identity": self.required_identity, "taxon": self.taxon,
            "trigger_type": self.trigger_type, "cap": self.cap, "bound": self.bound,
            "reference_accessions": list(self.reference_accessions),
            "reference_genotypes": list(self.reference_genotypes),
            "co_present_bacteria": list(self.co_present_bacteria),
            "required_context": list(self.required_context), "extrapolations": list(self.extrapolations),
        }


@dataclass(frozen=True, slots=True)
class Registry:
    taxa: tuple[TaxonEntry, ...]
    rules: tuple[Rule, ...]
    _index: Mapping[str, TaxonEntry] = field(default_factory=dict, repr=False)

    def lookup(self, name: str) -> TaxonEntry | None:
        """By accepted name or synonym, case-insensitive, underscores allowed."""
        key = _norm(name)
        hit = self._index.get(key)
        if hit is not None:
            return hit
        # Genus fallback for genus-level entries (Aspergillus, Penicillium, Fusarium).
        genus = key.split(" ")[0]
        entry = self._index.get(genus)
        return entry if entry is not None and entry.rank == "genus" else None

    @property
    def policy_lock_id(self) -> str:
        return policy_lock([r.to_json() for r in self.rules])

    def rule(self, rule_id: str) -> Rule | None:
        return next((r for r in self.rules if r.rule_id == rule_id), None)


def _norm(name: str) -> str:
    return re.sub(r"\s+", " ", name.replace("_", " ").strip()).lower()


@lru_cache(maxsize=1)
def load(path: Path = REGISTRY) -> Registry:
    data = yaml.safe_load(Path(path).read_text()) or {}
    taxa: list[TaxonEntry] = []
    for t in data.get("taxa") or []:
        claims = tuple(
            Claim(
                claim_id=str(c["claim_id"]), label=str(c.get("label") or "insufficient_evidence"),
                direction=str(c.get("direction") or "context"), endpoint=str(c.get("endpoint") or ""),
                host=str(c.get("host") or ""), site=str(c.get("site") or ""), design=str(c.get("design") or ""),
                source=str(c.get("source") or ""), limit=" ".join(str(c.get("limit") or "").split()),
                detail=" ".join(str(c.get("detail") or "").split()),
            )
            for c in (t.get("evidence") or [])
        )
        taxa.append(TaxonEntry(
            name=str(t["name"]), taxid=t.get("taxid"), rank=str(t.get("rank") or "species"),
            synonyms=tuple(str(s) for s in (t.get("synonyms") or [])), category=str(t.get("category") or "unknown"),
            origin_default=str(t.get("origin_default") or "unknown"), role=" ".join(str(t.get("role") or "").split()),
            evidence=claims, must_not_infer=" ".join(str(t.get("must_not_infer") or "").split()),
            complex=t.get("complex"),
        ))
    rules = tuple(
        Rule(
            rule_id=str(r["rule_id"]), direction=str(r["direction"]), biological_claim_id=str(r["biological_claim_id"]),
            source_ids=tuple(str(s) for s in (r.get("source_ids") or [])),
            required_identity=str(r.get("required_identity") or "species"), taxon=r.get("taxon"),
            trigger_type=str(r.get("trigger_type") or "presence"), cap=float(r.get("cap") or 0.0),
            bound=bool(r.get("bound", True)),
            reference_accessions=tuple(str(a) for a in (r.get("reference_accessions") or [])),
            reference_genotypes=tuple(str(g) for g in (r.get("reference_genotypes") or [])),
            co_present_bacteria=tuple(str(b) for b in (r.get("co_present_bacteria") or [])),
            required_context=tuple(str(c) for c in (r.get("required_context") or [])),
            extrapolations=tuple(str(e) for e in (r.get("extrapolations") or [])), note=" ".join(str(r.get("note") or "").split()),
        )
        for r in (data.get("rules") or [])
    )
    index: dict[str, TaxonEntry] = {}
    for t in taxa:
        index[_norm(t.name)] = t
        for s in t.synonyms:
            index.setdefault(_norm(s), t)
    return Registry(taxa=tuple(taxa), rules=rules, _index=index)


# --------------------------------------------------------------------------- #
# rule evaluation
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class SupportedTaxon:
    """What the evaluator needs to know about one supported fungal finding."""

    finding_id: str
    name: str                       # accepted species (or complex) name
    rank: str                       # species | complex | genus
    detection_state: str            # supported | provisional | ambiguous_complex
    strain_resolution: str          # spec §13 resolution values
    reference_accessions: tuple[str, ...] = ()
    reference_genotypes: tuple[str, ...] = ()
    #: Reference genotypes / accessions still compatible with the observed
    #: discriminating evidence (the equivalence set, or every panel member
    #: when nothing was callable).
    compatible_genotypes: tuple[str, ...] = ()
    compatible_accessions: tuple[str, ...] = ()
    strain_analysis_status: str = "not_assessed"


def evaluate(
    registry: Registry,
    supported: Sequence[SupportedTaxon],
    *,
    bacteria_present: Iterable[str] = (),
    context: Mapping[str, Any] | None = None,
) -> list[Contribution]:
    """Every rule in the manifest, evaluated against one sample.

    Returns one Contribution per rule, active or not, so the trace shows why
    each rule did or did not fire. Unbound rules are returned inactive with
    the reason from the manifest.
    """
    bact = {_norm(b) for b in bacteria_present}
    ctx = dict(context or {})
    by_name: dict[str, list[SupportedTaxon]] = {}
    for s in supported:
        if s.detection_state not in ("supported", "provisional", "ambiguous_complex"):
            continue
        entry = registry.lookup(s.name)
        key = _norm(entry.name) if entry is not None else _norm(s.name)
        by_name.setdefault(key, []).append(s)

    out: list[Contribution] = []
    for r in registry.rules:
        base = {
            "rule_id": r.rule_id, "direction": r.direction, "biological_claim_id": r.biological_claim_id,
            "source_ids": r.source_ids, "required_identity": r.required_identity, "trigger_type": r.trigger_type,
            "base_cap": r.cap, "required_features": (), "required_context": r.required_context,
            "extrapolations": r.extrapolations,
        }
        if not r.bound:
            out.append(Contribution(**base, activation=0.0, alternatives=(0.0,), activation_status="inactive",
                                    reason=f"unbound rule: {r.note}"))
            continue
        hits = [s for s in by_name.get(_norm(r.taxon or ""), []) if s.detection_state == "supported"]
        if not hits:
            out.append(Contribution(**base, activation=0.0, alternatives=(0.0,), activation_status="inactive",
                                    reason=f"{r.taxon} not supported in this sample"))
            continue
        fids = tuple(h.finding_id for h in hits)

        if r.trigger_type == "presence":
            if any(h.rank != "species" for h in hits) and all(h.rank != "species" for h in hits):
                out.append(Contribution(**base, activation=0.0, alternatives=(0.0,), activation_status="inactive",
                                        finding_ids=fids, reason="only a complex or higher-rank finding; species not identified"))
                continue
            out.append(Contribution(**base, activation=1.0, alternatives=(1.0,), activation_status="active",
                                    finding_ids=fids, reason=f"{r.taxon} supported at species level"))
            continue

        if r.trigger_type == "co_presence":
            missing = [b for b in r.co_present_bacteria if _norm(b) not in bact]
            if missing:
                out.append(Contribution(**base, activation=0.0, alternatives=(0.0,), activation_status="inactive",
                                        finding_ids=fids, reason=f"co-present bacteria not all detected: {', '.join(missing)}"))
            else:
                out.append(Contribution(**base, activation=1.0, alternatives=(1.0,), activation_status="active",
                                        finding_ids=fids, reason=f"{r.taxon} with {', '.join(r.co_present_bacteria)} all detected"))
            continue

        if r.trigger_type == "genotype":
            targets = set(r.reference_accessions) | set(r.reference_genotypes)
            resolved = [h for h in hits if h.strain_resolution == "reference_genotype"
                        and (set(h.reference_accessions) | set(h.reference_genotypes)) & targets]
            if resolved:
                if r.required_context and any(not ctx.get(c) for c in r.required_context):
                    out.append(Contribution(**base, activation=0.0, alternatives=(0.0, 1.0), activation_status="unresolved",
                                            finding_ids=fids, reason="genotype resolved but required context is not documented"))
                else:
                    out.append(Contribution(**base, activation=1.0, alternatives=(1.0,), activation_status="active",
                                            finding_ids=fids, reason="reference genotype resolved against near neighbours"))
                continue
            # Not resolved to the studied genotype. Is it still compatible?
            compatible = [h for h in hits if (set(h.compatible_accessions) | set(h.compatible_genotypes)) & targets]
            if r.cap <= 0:
                reason = "rule carries no points by policy; the finding stays visible"
            elif compatible:
                reason = ("studied genotype remains compatible with the observed discriminating evidence; "
                          "enters the sensitivity range only")
            elif any(h.strain_analysis_status == "complete" for h in hits):
                reason = "studied genotype excluded by discriminating loci"
            else:
                reason = "strain analysis did not complete; studied genotype not assessed"
            if compatible and r.cap > 0:
                out.append(Contribution(**base, activation=0.0, alternatives=(0.0, 1.0), activation_status="unresolved",
                                        finding_ids=fids, reason=reason))
            else:
                out.append(Contribution(**base, activation=0.0, alternatives=(0.0,), activation_status="inactive",
                                        finding_ids=fids, reason=reason))
            continue

        out.append(Contribution(**base, activation=0.0, alternatives=(0.0,), activation_status="inactive",
                                finding_ids=fids, reason=f"trigger type {r.trigger_type} needs a qualified reference"))
    return out


__all__ = [
    "CATEGORY_LABELS", "EVIDENCE_LABELS", "Claim", "Registry", "Rule", "SupportedTaxon", "TaxonEntry",
    "evaluate", "load",
]
