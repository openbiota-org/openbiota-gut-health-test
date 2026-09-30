"""What is known to act against a named organism, and what it is.

A pathogen section that names an organism and stops there leaves the reader
with the worst half of the information: enough to be frightened, not enough
to do anything. This module carries the other half — what the organism is,
whether carrying it is ordinary, what would make the amount concerning, and
which non-prescription agents have credible evidence of activity against it.

Three rules shape the content, and they are enforced on load rather than
trusted:

1. **Evidence tier is never optional.** An in-vitro finding and a randomised
   trial are both worth knowing and are not the same claim. Every agent
   carries its tier, and the report prints it next to the agent's name.
2. **The substance has to reach the gut.** An agent destroyed by stomach acid
   or absorbed before the colon cannot act on a colonic organism, however
   good the culture data. Each entry states why it arrives intact.
3. **An empty list is a real answer.** For an ordinary gut resident the
   honest recommendation is to leave it alone, and `no_agents_reason` says
   so. Padding those entries with plausible-sounding protocols would be the
   most harmful thing this module could do.

Nothing here is medical advice, and the loader refuses content that omits
the caution or the citation that would let a reader check it.
"""

from __future__ import annotations

import hashlib
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

import yaml

from openbiota.pathogens.schema import PathogenError

__all__ = [
    "AgentAdvice",
    "AgentLibrary",
    "EVIDENCE_TIERS",
    "OrganismAdvice",
    "TIER_WORDS",
    "load_agent_library",
]

#: Agent categories. Frozen: an unknown kind is a typo, not a new category.
AGENT_KINDS: Final[frozenset[str]] = frozenset(
    {"probiotic", "botanical", "nutrient", "food"}
)

#: Evidence tiers, strongest first. Order is meaningful: the report sorts by
#: it, so a human trial is never listed below a culture experiment.
EVIDENCE_TIERS: Final[tuple[str, ...]] = (
    "human_rct",
    "human_observational",
    "animal",
    "in_vitro",
)

#: How each tier is described to a reader. The wording is deliberately plain
#: about how far a culture result generalises, because that is exactly where
#: supplement claims usually overreach.
TIER_WORDS: Final[Mapping[str, str]] = {
    "human_rct": "randomised human trial",
    "human_observational": "human study, not randomised",
    "animal": "animal study only",
    "in_vitro": "laboratory culture only — no human evidence",
}

_REQUIRED_AGENT_FIELDS: Final[tuple[str, ...]] = (
    "name", "kind", "evidence_tier", "effect", "dose_studied",
    "reaches_gut", "citation_title", "citation_url", "caution",
)

_REQUIRED_ORGANISM_FIELDS: Final[tuple[str, ...]] = (
    "display_name", "what_it_is", "why_it_matters", "normal_carriage",
    "overgrowth_signal",
)


@dataclass(frozen=True, slots=True)
class AgentAdvice:
    """One agent with credible evidence against one organism."""

    name: str
    kind: str
    evidence_tier: str
    effect: str
    dose_studied: str
    reaches_gut: str
    citation_title: str
    citation_url: str
    caution: str

    @property
    def tier_rank(self) -> int:
        return EVIDENCE_TIERS.index(self.evidence_tier)

    @property
    def tier_words(self) -> str:
        return TIER_WORDS[self.evidence_tier]

    @property
    def is_human_evidence(self) -> bool:
        return self.evidence_tier.startswith("human")

    def to_json(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "kind": self.kind,
            "evidence_tier": self.evidence_tier,
            "effect": self.effect,
            "dose_studied": self.dose_studied,
            "reaches_gut": self.reaches_gut,
            "citation_title": self.citation_title,
            "citation_url": self.citation_url,
            "caution": self.caution,
        }


@dataclass(frozen=True, slots=True)
class OrganismAdvice:
    """Everything the report can say about one organism beyond its name."""

    seed_id: str
    display_name: str
    what_it_is: str
    why_it_matters: str
    normal_carriage: str
    overgrowth_signal: str
    agents: tuple[AgentAdvice, ...] = ()
    no_agents_reason: str = ""
    lane: str = ""

    @property
    def has_agents(self) -> bool:
        return bool(self.agents)

    @property
    def best_tier(self) -> str | None:
        return self.agents[0].evidence_tier if self.agents else None

    @property
    def human_backed(self) -> tuple[AgentAdvice, ...]:
        return tuple(a for a in self.agents if a.is_human_evidence)

    def to_json(self) -> dict[str, Any]:
        return {
            "seed_id": self.seed_id,
            "display_name": self.display_name,
            "what_it_is": self.what_it_is,
            "why_it_matters": self.why_it_matters,
            "normal_carriage": self.normal_carriage,
            "overgrowth_signal": self.overgrowth_signal,
            "agents": [a.to_json() for a in self.agents],
            "no_agents_reason": self.no_agents_reason,
            "lane": self.lane,
        }


@dataclass(frozen=True, slots=True)
class AgentLibrary:
    """Advice keyed by catalogue seed id."""

    by_seed: Mapping[str, OrganismAdvice] = field(default_factory=dict)
    lanes: tuple[str, ...] = ()

    def get(self, seed_id: str) -> OrganismAdvice | None:
        return self.by_seed.get(seed_id)

    @property
    def n_agents(self) -> int:
        return sum(len(o.agents) for o in self.by_seed.values())

    def __len__(self) -> int:
        return len(self.by_seed)

    def __bool__(self) -> bool:
        return bool(self.by_seed)


def _clean(text: Any) -> str:
    """YAML block scalars arrive with newlines; the report wants a sentence."""
    return " ".join(str(text or "").split())


def _agent(seed_id: str, raw: Mapping[str, Any]) -> AgentAdvice:
    for key in _REQUIRED_AGENT_FIELDS:
        if not str(raw.get(key) or "").strip():
            raise PathogenError(
                f"{seed_id}: agent {raw.get('name', '?')!r} is missing {key}. "
                "Every agent must carry its evidence tier, why it reaches the "
                "gut, a checkable citation and its caution."
            )
    kind = str(raw["kind"]).strip()
    tier = str(raw["evidence_tier"]).strip()
    if kind not in AGENT_KINDS:
        raise PathogenError(f"{seed_id}: unknown agent kind {kind!r}")
    if tier not in EVIDENCE_TIERS:
        raise PathogenError(f"{seed_id}: unknown evidence tier {tier!r}")
    url = str(raw["citation_url"]).strip()
    if not url.startswith("http"):
        raise PathogenError(f"{seed_id}: citation_url is not a URL: {url!r}")
    return AgentAdvice(
        name=_clean(raw["name"]),
        kind=kind,
        evidence_tier=tier,
        effect=_clean(raw["effect"]),
        dose_studied=_clean(raw["dose_studied"]),
        reaches_gut=_clean(raw["reaches_gut"]),
        citation_title=_clean(raw["citation_title"]),
        citation_url=url,
        caution=_clean(raw["caution"]),
    )


def _organism(seed_id: str, raw: Mapping[str, Any], lane: str) -> OrganismAdvice:
    for key in _REQUIRED_ORGANISM_FIELDS:
        if not str(raw.get(key) or "").strip():
            raise PathogenError(f"{seed_id}: missing {key}")
    agents = [_agent(seed_id, a) for a in (raw.get("agents") or ())]
    # Strongest evidence first, then alphabetically, so the ordering is
    # deterministic and a culture result never outranks a trial.
    agents.sort(key=lambda a: (a.tier_rank, a.name.lower()))
    reason = _clean(raw.get("no_agents_reason"))
    if not agents and not reason:
        raise PathogenError(
            f"{seed_id}: no agents and no `no_agents_reason`. An empty list is "
            "a legitimate answer for an ordinary resident, but it has to say why."
        )
    return OrganismAdvice(
        seed_id=seed_id,
        display_name=_clean(raw["display_name"]),
        what_it_is=_clean(raw["what_it_is"]),
        why_it_matters=_clean(raw["why_it_matters"]),
        normal_carriage=_clean(raw["normal_carriage"]),
        overgrowth_signal=_clean(raw["overgrowth_signal"]),
        agents=tuple(agents),
        no_agents_reason=reason,
        lane=lane,
    )


def load_agent_library(
    catalog_dir: Path,
    *,
    known_seed_ids: Sequence[str] | None = None,
) -> AgentLibrary:
    """Load every `agents_*.yaml` beside the catalogue.

    Absence is not an error: the section renders without the advice column
    when no file is installed. Malformed content *is* an error, because
    silently dropping half an entry would leave the report stating an
    organism is treatable without saying with what.
    """
    from openbiota.fastcache import cached_load

    directory = Path(catalog_dir)
    files = sorted(directory.glob("agents_*.yaml"))
    if not files:
        return AgentLibrary()
    # The allow-list is part of what the result depends on; a different list
    # is a different library, so it is folded into the cache name.
    allowed_key = "" if known_seed_ids is None else hashlib.sha256(
        "\n".join(sorted(known_seed_ids)).encode()
    ).hexdigest()[:12]
    return cached_load(
        f"agents{'-' + allowed_key if allowed_key else ''}", files,
        lambda: _load_agent_library_uncached(files, known_seed_ids),
    )


def _load_agent_library_uncached(
    files: Sequence[Path], known_seed_ids: Sequence[str] | None
) -> AgentLibrary:
    allowed = set(known_seed_ids) if known_seed_ids is not None else None
    by_seed: dict[str, OrganismAdvice] = {}
    lanes: list[str] = []
    for path in files:
        try:
            doc = yaml.safe_load(path.read_text()) or {}
        except yaml.YAMLError as exc:
            raise PathogenError(f"{path.name}: {exc}") from exc
        lane = str(doc.get("lane") or path.stem)
        lanes.append(lane)
        organisms = doc.get("organisms") or {}
        if not isinstance(organisms, Mapping):
            raise PathogenError(f"{path.name}: `organisms` must be a mapping")
        for seed_id, raw in organisms.items():
            key = str(seed_id)
            if allowed is not None and key not in allowed:
                raise PathogenError(
                    f"{path.name}: {key} is not a catalogue seed id. Advice "
                    "keyed to a name no target carries would never be shown."
                )
            if key in by_seed:
                raise PathogenError(
                    f"{key} appears in two agent files ({by_seed[key].lane} "
                    f"and {lane}); one organism, one entry."
                )
            by_seed[key] = _organism(key, raw, lane)
    return AgentLibrary(by_seed=by_seed, lanes=tuple(lanes))


def self_test() -> int:
    """Checks that do not need the installed files."""
    checks = 0

    good = {
        "display_name": "Test organism",
        "what_it_is": "A test.",
        "why_it_matters": "It is a test.",
        "normal_carriage": "Never.",
        "overgrowth_signal": "None.",
        "agents": [
            {
                "name": "Culture agent", "kind": "botanical",
                "evidence_tier": "in_vitro", "effect": "Inhibited in culture.",
                "dose_studied": "not established", "reaches_gut": "Poorly absorbed.",
                "citation_title": "A paper", "citation_url": "https://example.org/a",
                "caution": "Not a treatment.",
            },
            {
                "name": "Trial agent", "kind": "probiotic",
                "evidence_tier": "human_rct", "effect": "Reduced recurrence.",
                "dose_studied": "500 mg twice daily", "reaches_gut": "Acid stable.",
                "citation_title": "A trial", "citation_url": "https://example.org/b",
                "caution": "Avoid if immunosuppressed.",
            },
        ],
    }
    org = _organism("x.y", good, "test")
    # Human evidence sorts above culture, whatever order it arrived in.
    assert org.agents[0].name == "Trial agent", "tiers must sort strongest first"
    assert org.best_tier == "human_rct"
    assert len(org.human_backed) == 1
    checks += 1

    # An empty list needs a stated reason; a missing field is refused.
    for broken, why in (
        ({**good, "agents": []}, "empty without reason"),
        ({**good, "agents": [{**good["agents"][0], "citation_url": "nope"}]},
         "citation not a url"),
        ({**good, "agents": [{**good["agents"][0], "evidence_tier": "vibes"}]},
         "unknown tier"),
        ({**good, "agents": [{**good["agents"][0], "caution": ""}]},
         "missing caution"),
        ({**good, "why_it_matters": ""}, "missing organism field"),
    ):
        try:
            _organism("x.y", broken, "test")
        except PathogenError:
            checks += 1
        else:  # pragma: no cover - guard
            raise AssertionError(f"accepted {why}")

    # An empty list *with* a reason is fine, and is not a failure state.
    left_alone = _organism(
        "x.z", {**good, "agents": [], "no_agents_reason": "An ordinary resident."},
        "test",
    )
    assert not left_alone.has_agents and left_alone.no_agents_reason
    assert left_alone.best_tier is None
    checks += 1

    # Every tier has reader-facing words, and the culture tier says so.
    assert set(TIER_WORDS) == set(EVIDENCE_TIERS)
    assert "no human evidence" in TIER_WORDS["in_vitro"]
    checks += 1

    # A missing directory is silence, not a crash.
    empty = load_agent_library(Path("/nonexistent-agents-dir"))
    assert not empty and empty.n_agents == 0
    checks += 1

    return checks


if __name__ == "__main__":  # pragma: no cover
    print(f"{self_test()} agent-library checks passed")
