"""Gut community type — A09, BUILD_SPEC_v0.8.3 §6.2.

Two things live here, and the separation is the point.

**The composition view is always available.** A genus-level composition over
non-overlapping partitions, normalised to one, with an explicit
`unresolved` share. It needs no model, no cohort and no clustering, and it
is what a reader actually looks at.

**The classifier is a frozen artifact or it is nothing.** Assignment is to
the nearest medoid of a model trained once, on a fixed cohort, and written
to disk with its evaluation report. Nothing here clusters at report time,
and nothing ranks a sample against whatever other samples happen to be in
the run. If the frozen model says its own clustering was unstable, this
module returns the composition and the distances and declines to name a
type - which the specification requires and which is the honest outcome
when the data do not support discrete groups.

Jensen-Shannon divergence with natural logs; PAM on `sqrt(JSD)`, which is a
metric where JSD itself is not. Soft similarities, where shown, are
`softmax(-d / tau)` with `tau` chosen on held-out reference data, and they
are similarities to reference groups - never probabilities of a disease.
"""

from __future__ import annotations

import json
import math
import re
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

MODEL_FILE: Final = Path(__file__).resolve().parents[2] / "extension" / "community_type_model.json"
METHOD_COMPOSITION: Final = "ext083.genus_composition/1.0"
METHOD_COMMUNITY_TYPE: Final = "ext083.community_type/1.0"

#: The label for everything the genus universe does not name. It is a
#: partition of the composition, not a residual to be dropped: a sample that
#: is 40% organisms the model never saw is a different sample from one that
#: is 40% Bacteroides, and the number has to be on the page.
UNRESOLVED: Final = "unresolved"

#: How close two distances have to be before an assignment is called mixed
#: rather than given a name, as a fraction of the nearest distance.
MIXED_MARGIN: Final = 0.05

_GTDB_SUFFIX: Final = re.compile(r"_[A-Z]{1,2}$")

#: Genus reclassifications, current name -> the name the reference cohort
#: uses. The frozen cohort was profiled with MetaPhlAn 3, whose taxonomy
#: predates a decade of transfers; without this table this sample's single
#: largest genus, *Phocaeicola* at 21%, is simply absent from the reference
#: and the distance is computed on two fifths of the sample.
#:
#: Only genuine transfers are listed: the type species moved from the old
#: genus to the new one, so the two names denote the same organisms. A
#: merely *related* genus is not a synonym and is not here - *Pseudo-
#: flavonifractor* stays out of *Flavonifractor*, and an unnamed MetaPhlAn
#: bin (`GGB...`) stays unresolved rather than being pushed into a
#: convenient neighbour.
GENUS_SYNONYMS: Final[Mapping[str, str]] = {
    # García-López et al. 2019, Front. Microbiol. 10:2083 - Bacteroidaceae split.
    "Phocaeicola": "Bacteroides",
    # Togo et al. 2018, New Microbes New Infect. - Ruminococcaceae/Lachnospiraceae transfer.
    "Mediterraneibacter": "Ruminococcus",
    # Rosero et al. 2016, IJSEM 66:768 - Eubacterium rectale -> Agathobacter rectalis.
    "Agathobacter": "Eubacterium",
    # Shetty et al. 2018, IJSEM 68:3741 - Eubacterium hallii -> Anaerobutyricum hallii.
    "Anaerobutyricum": "Eubacterium",
    # Haas & Blanchard 2020, IJSEM 70:23 - Clostridium clusters XIVa/XVIII split.
    "Enterocloster": "Clostridium",
    "Hungatella": "Clostridium",
    "Thomasclavelia": "Clostridium",
    "Erysipelatoclostridium": "Clostridium",
    "Tyzzerella": "Clostridium",
    "Clostridioides": "Clostridium",
    # Ueki et al. 2017, IJSEM 67:4146.
    "Anaerotignum": "Clostridium",
    # Zheng et al. 2020, IJSEM 70:2782 - the Lactobacillus genus split into 25.
    "Lactiplantibacillus": "Lactobacillus",
    "Limosilactobacillus": "Lactobacillus",
    "Lacticaseibacillus": "Lactobacillus",
    "Ligilactobacillus": "Lactobacillus",
    "Latilactobacillus": "Lactobacillus",
    "Lentilactobacillus": "Lactobacillus",
    "Companilactobacillus": "Lactobacillus",
    "Loigolactobacillus": "Lactobacillus",
    "Furfurilactobacillus": "Lactobacillus",
    "Schleiferilactobacillus": "Lactobacillus",
    "Levilactobacillus": "Lactobacillus",
    # Hitch et al. 2024, IJSEM - Prevotella copri -> Segatella copri. This one
    # decides an enterotype, so leaving it out would change the answer.
    "Segatella": "Prevotella",
    "Hoylesella": "Prevotella",
    "Leyella": "Prevotella",
    "Xylanibacter": "Prevotella",
    # Liu et al. 2021 / GTDB - Coprococcus transfers.
    "Allocoprococcus": "Coprococcus",
    # Ludwig et al. 2015 - Faecalibacterium / Gemmiger neighbourhood, type
    # species transferred out of Ruminococcus.
    "Agathobaculum": "Agathobaculum",
}


def canonical_genus(genus: str) -> str:
    """The genus under the reference cohort's taxonomy."""
    return GENUS_SYNONYMS.get(genus, genus)


class CommunityTypeError(RuntimeError):
    """A malformed frozen model, or a request that needs one and has none."""


# --------------------------------------------------------------------------- #
# composition — always available
# --------------------------------------------------------------------------- #


def genus_of(species: str) -> str:
    """The genus a species label belongs to, GTDB split suffix removed.

    `Prevotella_A copri` and `Prevotella copri` are one genus for the purpose
    of a composition chart; keeping them apart would split a single bar in
    two for a reason no reader could see.
    """
    text = str(species).replace("_", " ").strip()
    if not text:
        return ""
    first = text.split()[0]
    first = first.strip("[]")
    return _GTDB_SUFFIX.sub("", first)


@dataclass(frozen=True)
class Composition:
    """A genus composition over non-overlapping partitions, summing to one."""

    shares: dict[str, float]
    n_taxa: int
    n_genera: int
    total_percent: float
    denominator: str

    @property
    def unresolved_share(self) -> float:
        return self.shares.get(UNRESOLVED, 0.0)

    def top(self, n: int = 12) -> list[tuple[str, float]]:
        named = [(k, v) for k, v in self.shares.items() if k != UNRESOLVED]
        named.sort(key=lambda kv: (-kv[1], kv[0]))
        return named[:n]

    def to_json(self) -> dict[str, Any]:
        return {
            "method_id": METHOD_COMPOSITION,
            "shares": dict(sorted(self.shares.items(), key=lambda kv: (-kv[1], kv[0]))),
            "top_genera": [{"genus": g, "share": s} for g, s in self.top()],
            "unresolved_share": self.unresolved_share,
            "n_taxa": self.n_taxa,
            "n_genera": self.n_genera,
            "total_percent_of_sample": self.total_percent,
            "denominator": self.denominator,
            "partition_note": (
                "Non-overlapping partitions of one lane's abundance, normalised to one. "
                "`unresolved` is everything the genus universe does not name; it is a share "
                "of the same total, not a residual that has been dropped."
            ),
        }


def composition(
    abundances: Mapping[str, float],
    *,
    universe: Sequence[str] | None = None,
    denominator: str = "bacterial abundance in the primary lane",
) -> Composition:
    """Genus composition, with everything outside `universe` pooled as unresolved."""
    known = set(universe) if universe else None
    by_genus: dict[str, float] = {}
    total = 0.0
    n_taxa = 0
    for species, percent in abundances.items():
        value = float(percent or 0.0)
        if value <= 0:
            continue
        n_taxa += 1
        total += value
        genus = genus_of(species)
        key = genus if (known is None or genus in known) and genus else UNRESOLVED
        by_genus[key] = by_genus.get(key, 0.0) + value
    if total <= 0:
        return Composition({}, 0, 0, 0.0, denominator)
    shares = {k: v / total for k, v in by_genus.items()}
    return Composition(
        shares=shares, n_taxa=n_taxa,
        n_genera=sum(1 for k in shares if k != UNRESOLVED),
        total_percent=total, denominator=denominator,
    )


# --------------------------------------------------------------------------- #
# divergence
# --------------------------------------------------------------------------- #


def jensen_shannon(p: Sequence[float], q: Sequence[float]) -> float:
    """Jensen-Shannon divergence, natural logs, in nats.

    Zero where a component is zero in both; the usual `0 log 0 = 0`. The
    result is bounded by `ln 2`, so `sqrt(JSD)` is bounded by `sqrt(ln 2)`.
    """
    if len(p) != len(q):
        raise CommunityTypeError(f"vectors of different length: {len(p)} and {len(q)}")
    total = 0.0
    for a, b in zip(p, q, strict=True):
        m = 0.5 * (a + b)
        if m <= 0:
            continue
        if a > 0:
            total += 0.5 * a * math.log(a / m)
        if b > 0:
            total += 0.5 * b * math.log(b / m)
    return max(0.0, total)


def root_jsd(p: Sequence[float], q: Sequence[float]) -> float:
    """`sqrt(JSD)`, which satisfies the triangle inequality where JSD does not."""
    return math.sqrt(jensen_shannon(p, q))


def under_reference_taxonomy(shares: Mapping[str, float]) -> dict[str, float]:
    """A composition restated in the reference cohort's genus names."""
    out: dict[str, float] = {}
    for genus, share in shares.items():
        out[canonical_genus(genus)] = out.get(canonical_genus(genus), 0.0) + float(share)
    return out


def vectorise(
    shares: Mapping[str, float], features: Sequence[str], *, apply_synonyms: bool = True
) -> list[float]:
    """A composition as a vector over the model's fixed features, renormalised.

    Renormalisation matters: a sample whose named genera cover 60% of it must
    be compared with the reference on the same basis the reference used, or
    every distance is inflated by the coverage difference rather than by the
    composition.
    """
    source = under_reference_taxonomy(shares) if apply_synonyms else shares
    raw = [max(0.0, float(source.get(name, 0.0))) for name in features]
    total = sum(raw)
    if total <= 0:
        return [0.0] * len(features)
    return [v / total for v in raw]


# --------------------------------------------------------------------------- #
# the frozen model
# --------------------------------------------------------------------------- #


@dataclass(frozen=True)
class CommunityTypeModel:
    """A model trained once, evaluated once, and written to disk."""

    model_id: str
    features: tuple[str, ...]
    medoids: tuple[tuple[float, ...], ...]
    labels: tuple[str, ...]
    distinguishing_taxa: tuple[tuple[str, ...], ...]
    k: int
    tau: float
    stable: bool
    evaluation: Mapping[str, Any]
    cohort: Mapping[str, Any]
    reference_coverage: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def load(cls, path: Path | None = None) -> CommunityTypeModel:
        file = path or MODEL_FILE
        if not file.is_file():
            raise CommunityTypeError(f"no frozen community-type model at {file}")
        try:
            raw = json.loads(file.read_text(encoding="utf-8"))
            model = cls(
                model_id=str(raw["model_id"]),
                features=tuple(str(f) for f in raw["features"]),
                medoids=tuple(tuple(float(x) for x in row) for row in raw["medoids"]),
                labels=tuple(str(x) for x in raw["labels"]),
                distinguishing_taxa=tuple(
                    tuple(str(t) for t in row) for row in raw["distinguishing_taxa"]
                ),
                k=int(raw["k"]), tau=float(raw["tau"]), stable=bool(raw["stable"]),
                evaluation=raw.get("evaluation") or {}, cohort=raw.get("cohort") or {},
                reference_coverage=raw.get("reference_feature_coverage") or {},
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise CommunityTypeError(f"malformed community-type model: {exc}") from exc
        if len(model.medoids) != model.k or len(model.labels) != model.k:
            raise CommunityTypeError(
                f"model declares k={model.k} but carries {len(model.medoids)} medoids "
                f"and {len(model.labels)} labels"
            )
        for row in model.medoids:
            if len(row) != len(model.features):
                raise CommunityTypeError("a medoid does not span the feature set")
        return model


@dataclass(frozen=True)
class Assignment:
    """Where a sample sits relative to every reference group."""

    model_id: str
    distances: tuple[tuple[str, float], ...]
    similarities: tuple[tuple[str, float], ...]
    nearest: str | None
    margin: float | None
    state: str
    coverage: float
    reference_coverage: Mapping[str, Any] = field(default_factory=dict)
    #: (current name, reference name, share) for every genus restated to
    #: reach the reference's taxonomy. Shown, not hidden: a reader comparing
    #: this page with the organism list must be able to see why a genus
    #: appears under an older name here.
    renamed_genera: tuple[tuple[str, str, float], ...] = ()
    distinguishing_taxa: tuple[str, ...] = ()
    note: str = ""

    def to_json(self) -> dict[str, Any]:
        return {
            "method_id": METHOD_COMMUNITY_TYPE,
            "model_id": self.model_id,
            "state": self.state,
            "nearest_group": self.nearest,
            "margin_to_next": self.margin,
            "distances": [{"group": g, "root_jsd": d} for g, d in self.distances],
            "similarities": [{"group": g, "similarity": s} for g, s in self.similarities],
            "distinguishing_taxa": list(self.distinguishing_taxa),
            "feature_coverage": self.coverage,
            "reference_feature_coverage": dict(self.reference_coverage),
            "renamed_for_comparison": [
                {"as_reported": a, "as_in_reference": b, "share": c}
                for a, b, c in self.renamed_genera
            ],
            "coverage_note": (
                "the share of this sample the classifier's genus universe names. Distances are "
                "computed on that share, renormalised, exactly as the reference was, so a "
                "coverage near the reference range does not weaken the comparison"
            ),
            "note": self.note,
            "limitations": [
                "A resemblance to a reference group, not a diagnosis and not a disease "
                "probability.",
                "A Bacteroides-dominant result does not establish a high-protein or high-fat "
                "intake.",
                "A relative composition resembling a published low-load group does not "
                "establish low microbial biomass; that needs an absolute-load assay.",
            ],
        }


def assign(
    comp: Composition, model: CommunityTypeModel | None = None, *, path: Path | None = None
) -> Assignment:
    """Place a composition against the frozen model's medoids.

    Every distance is reported, not only the nearest. Where the model itself
    recorded unstable clustering, or where two groups are within
    `MIXED_MARGIN`, the sample is `mixed` and carries no group name.
    """
    model = model or CommunityTypeModel.load(path)
    restated = under_reference_taxonomy(comp.shares)
    covered = sum(restated.get(f, 0.0) for f in model.features)
    feature_set = set(model.features)
    renamed = tuple(sorted(
        (
            (genus, canonical_genus(genus), share)
            for genus, share in comp.shares.items()
            if canonical_genus(genus) != genus and canonical_genus(genus) in feature_set
        ),
        key=lambda row: -row[2],
    ))
    vector = vectorise(comp.shares, model.features)
    if sum(vector) <= 0:
        return Assignment(
            model_id=model.model_id, distances=(), similarities=(), nearest=None, margin=None,
            state="unavailable", coverage=covered, reference_coverage=model.reference_coverage,
            renamed_genera=renamed,
            note=("none of this sample's abundance falls in the model's genus universe, so no "
                  "distance can be computed"),
        )
    pairs = sorted(
        ((label, root_jsd(vector, medoid))
         for label, medoid in zip(model.labels, model.medoids, strict=True)),
        key=lambda kv: (kv[1], kv[0]),
    )
    exponents = [-d / model.tau for _, d in pairs]
    peak = max(exponents)
    weights = [math.exp(e - peak) for e in exponents]
    total = sum(weights) or 1.0
    similarities = tuple((label, w / total) for (label, _), w in zip(pairs, weights, strict=True))

    nearest_label, nearest_d = pairs[0]
    margin = (pairs[1][1] - nearest_d) / nearest_d if len(pairs) > 1 and nearest_d > 0 else None
    if not model.stable:
        return Assignment(
            model_id=model.model_id, distances=tuple(pairs), similarities=similarities,
            nearest=None, margin=margin, state="continuous_only", coverage=covered,
            reference_coverage=model.reference_coverage, renamed_genera=renamed,
            note=("the frozen model's clustering was not stable enough to name discrete types, "
                  "so the composition and its distances are reported instead of a category"),
        )
    if margin is not None and margin < MIXED_MARGIN:
        return Assignment(
            model_id=model.model_id, distances=tuple(pairs), similarities=similarities,
            nearest=None, margin=margin, state="mixed", coverage=covered,
            reference_coverage=model.reference_coverage, renamed_genera=renamed,
            note=(f"this sample sits between {pairs[0][0]} and {pairs[1][0]}; the two are within "
                  f"{margin:.1%} of each other, which is too close to name one"),
        )
    index = model.labels.index(nearest_label)
    return Assignment(
        model_id=model.model_id, distances=tuple(pairs), similarities=similarities,
        nearest=nearest_label, margin=margin, state="assigned", coverage=covered,
        reference_coverage=model.reference_coverage, renamed_genera=renamed,
        distinguishing_taxa=model.distinguishing_taxa[index],
        note=(f"closest to the {nearest_label} reference group, named for the genera that "
              "distinguish it in the reference cohort"),
    )


# --------------------------------------------------------------------------- #
# training — run once, offline, by `scripts/train_community_type.py`
# --------------------------------------------------------------------------- #


@dataclass
class TrainingReport:
    k_evaluated: list[int] = field(default_factory=list)
    silhouette: dict[int, float] = field(default_factory=dict)
    stability: dict[int, float] = field(default_factory=dict)
    selected_k: int | None = None
    stable: bool = False
    reason: str = ""

    def to_json(self) -> dict[str, Any]:
        return {
            "k_evaluated": self.k_evaluated,
            "held_out_silhouette": {str(k): v for k, v in self.silhouette.items()},
            "study_stratified_bootstrap_stability": {str(k): v for k, v in self.stability.items()},
            "selected_k": self.selected_k,
            "stable": self.stable,
            "selection_rule": (
                "highest held-out silhouette among k with bootstrap stability at or above the "
                "declared threshold; if no k clears it, the model is marked unstable and "
                "assignment returns composition and distances only"
            ),
            "reason": self.reason,
        }


def pam(
    distance: Sequence[Sequence[float]], k: int, *, max_iter: int = 60, seed: int = 0
) -> list[int]:
    """Partitioning around medoids. Deterministic: k-medoids++ from a fixed seed.

    Written out rather than imported so that the frozen artifact does not
    depend on a clustering library's version-to-version tie-breaking.
    """
    import random  # noqa: PLC0415

    n = len(distance)
    if k >= n:
        return list(range(n))
    rng = random.Random(seed)
    medoids = [rng.randrange(n)]
    while len(medoids) < k:
        # k-means++ style: the farthest-from-any-medoid point wins, in proportion.
        weights = [min(distance[i][m] for m in medoids) ** 2 for i in range(n)]
        total = sum(weights)
        if total <= 0:
            medoids.append(next(i for i in range(n) if i not in medoids))
            continue
        pick = rng.random() * total
        running = 0.0
        for i, w in enumerate(weights):
            running += w
            if running >= pick and i not in medoids:
                medoids.append(i)
                break
        else:
            medoids.append(next(i for i in range(n) if i not in medoids))
    medoids.sort()
    for _ in range(max_iter):
        labels = [min(range(k), key=lambda j: distance[i][medoids[j]]) for i in range(n)]
        moved = False
        for j in range(k):
            members = [i for i in range(n) if labels[i] == j]
            if not members:
                continue
            best = min(members, key=lambda c: sum(distance[c][i] for i in members))
            if best != medoids[j]:
                medoids[j] = best
                moved = True
        if not moved:
            break
    medoids.sort()
    return medoids


def silhouette(distance: Sequence[Sequence[float]], labels: Sequence[int], k: int) -> float:
    """Mean silhouette width. Singleton clusters contribute zero, not one."""
    n = len(distance)
    scores: list[float] = []
    for i in range(n):
        own = [j for j in range(n) if labels[j] == labels[i] and j != i]
        if not own:
            scores.append(0.0)
            continue
        a = sum(distance[i][j] for j in own) / len(own)
        others = []
        for c in range(k):
            members = [j for j in range(n) if labels[j] == c]
            if c == labels[i] or not members:
                continue
            others.append(sum(distance[i][j] for j in members) / len(members))
        if not others:
            scores.append(0.0)
            continue
        b = min(others)
        scores.append((b - a) / max(a, b) if max(a, b) > 0 else 0.0)
    return sum(scores) / len(scores) if scores else 0.0


def distinguishing(
    vectors: Sequence[Sequence[float]], labels: Sequence[int], features: Sequence[str],
    cluster: int, *, top: int = 3,
) -> tuple[str, ...]:
    """The genera that separate one cluster from the rest, by mean difference."""
    inside = [v for v, lab in zip(vectors, labels, strict=True) if lab == cluster]
    outside = [v for v, lab in zip(vectors, labels, strict=True) if lab != cluster]
    if not inside or not outside:
        return ()
    scored: list[tuple[float, str]] = []
    for idx, name in enumerate(features):
        mean_in = sum(v[idx] for v in inside) / len(inside)
        mean_out = sum(v[idx] for v in outside) / len(outside)
        scored.append((mean_in - mean_out, name))
    scored.sort(reverse=True)
    return tuple(name for diff, name in scored[:top] if diff > 0)


def genus_matrix(
    samples: Iterable[Mapping[str, float]], *, min_prevalence: float = 0.10,
    min_mean_share: float = 0.001,
) -> tuple[list[str], list[list[float]]]:
    """A fixed genus universe and the cohort's vectors over it.

    The universe is fixed here, once, and written into the artifact: a genus
    that enters or leaves it later would silently change every distance.
    """
    compositions = [composition(s).shares for s in samples]
    if not compositions:
        return [], []
    counts: dict[str, int] = {}
    totals: dict[str, float] = {}
    for shares in compositions:
        for genus, share in shares.items():
            if genus == UNRESOLVED or share <= 0:
                continue
            counts[genus] = counts.get(genus, 0) + 1
            totals[genus] = totals.get(genus, 0.0) + share
    n = len(compositions)
    features = sorted(
        g for g in counts
        if counts[g] / n >= min_prevalence and totals[g] / n >= min_mean_share
    )
    vectors = [vectorise(shares, features) for shares in compositions]
    return features, vectors


__all__ = [
    "METHOD_COMMUNITY_TYPE",
    "METHOD_COMPOSITION",
    "MIXED_MARGIN",
    "MODEL_FILE",
    "UNRESOLVED",
    "Assignment",
    "CommunityTypeError",
    "CommunityTypeModel",
    "Composition",
    "TrainingReport",
    "GENUS_SYNONYMS",
    "assign",
    "canonical_genus",
    "composition",
    "distinguishing",
    "genus_matrix",
    "genus_of",
    "jensen_shannon",
    "pam",
    "root_jsd",
    "silhouette",
    "under_reference_taxonomy",
    "vectorise",
]
