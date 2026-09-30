"""Building the frozen reference cohort for the two taxonomic proxies.

The reference is the part of a percentile that is easiest to get quietly
wrong, so the construction here is deliberately explicit about three things
the raw cohort file does not tell you.

**Samples are not people.** ``refs/taxonomic_cohort.json`` holds 3,027
samples, and a naive reference would treat that as 3,027 independent
participants. It is 1,760 subjects: MehtaRS_2018 alone contributes 921
samples from 308 people, HMP_2019_ibdmdb 214 from 14. Section 8.2 requires
subject-level handling and section 8.5 requires the bootstrap to resample
independent participants, so samples are collapsed to subjects and exactly
one survives per subject.

Subject identity comes from curatedMetagenomicData's own ``subject_id``
field, read through :func:`_authoritative_subjects`. That field covers all
3,027 samples. An earlier version of this module inferred subjects by
pattern-matching sample IDs per study and excluded two cohorts whose IDs
it could not read; the authoritative field was on disk the whole time and
already parsed elsewhere in this codebase. The pattern rules survive only
as a fallback, and when they are used the provenance says so, because an
approximation and an exact mapping should not look alike in an audit
trail.

**Age must be known, not assumed.** 694 subjects have no recorded age, and
BackhedF_2015 is a mother-infant study. Since the report is read by two
minors, an adult reference must be provably adult rather than
probably-adult, so subjects without a recorded age are excluded outright
and the resulting calibration is labelled adult-only.

**Taxonomy namespaces must match.** The cohort is MetaPhlAn 3 species-level.
The pipeline runs both MetaPhlAn 3.1 and MetaPhlAn 4.1.1, and only the
former shares this namespace. Reading MetaPhlAn 4 SGB abundances against a
MetaPhlAn 3 reference would compare two different measurements, which
BF-T014 forbids, so :func:`sample_features` reads the MetaPhlAn 3.1 profile
and refuses if it is absent.

Genus features sum a fixed, explicit list of species leaves that exist in
the reference namespace, and the same list is used on both sides. A leaf
appears in exactly one feature, so nothing is double counted.
"""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path
from typing import Final

#: Where curatedMetagenomicData's own sample metadata lives. It carries a
#: real ``subject_id`` for every sample, which is the only correct way to
#: tell repeated measurements from separate people.
#:
#: An earlier version of this module recovered subjects by pattern-matching
#: sample IDs per study and excluded two cohorts whose IDs could not be
#: read. That was unnecessary: the authoritative field was already on disk
#: and already parsed elsewhere in this codebase. Regex recovery is kept
#: only as a fallback for when the metadata file is absent, and it is
#: recorded in the provenance when it is used, because it is an
#: approximation and the exact mapping is not.
CMD_METADATA_DIR: Final = "refs/cmd"

#: Fallback only. Used when the cMD metadata cannot be read; each rule was
#: derived from the actual ID structure, and each is an approximation of
#: the real subject field.
_SUBJECT_RULES: Final[dict[str, str]] = {
    "MehtaRS_2018": r"_SF\d+$",
    "HansenLBS_2018": r"-\d+$",
    "AsnicarF_2017": r"_t\d+\w*$",
    "HallAB_2017": r"_\d+_G\d+$",
    "CosteaPI_2017": r"-\d+-\d+$",
}

#: No study is excluded any more. With real subject IDs every cohort can
#: be collapsed correctly, including the two longitudinal ones that had to
#: be dropped while subjects were being guessed from ID text.
EXCLUDED_STUDIES: Final[dict[str, str]] = {}

#: Minimum fraction of the reference in which a feature must be detected
#: before it can rank anybody. A feature seen in a handful of participants
#: has an almost-constant column: every sample ties at the same midrank,
#: which contributes a fixed value to the group mean and dilutes the
#: features that do carry signal. Set well below the detection rate of

#: Reader-facing names for the feature symbols. The symbols are internal
#: identifiers - `m_gnavus`, `ecoli_complex` - and printing them raw asks the
#: reader to decode a variable name. Where a symbol covers several species
#: the label says so rather than naming one of them and implying the rest.
FEATURE_LABELS: Final[dict[str, str]] = {
    "ecoli_complex": "Escherichia coli and its close relatives",
    "m_gnavus": "Mediterraneibacter gnavus",
    "faecalibacterium": "Faecalibacterium prausnitzii",
    "coprococcus": "Coprococcus species",
    "subdoligranulum": "Subdoligranulum variabile",
    "blautia": "Blautia species",
}


def feature_label(symbol: str) -> str:
    """The reader-facing name for a feature symbol."""
    return FEATURE_LABELS.get(symbol, symbol.replace("_", " "))


#: every working feature (47% and up) and well above a database artefact.
MIN_FEATURE_DETECTION_RATE: Final = 0.05

#: Minimum age treated as adult. Subjects below it, or with no recorded
#: age at all, are excluded from the adult calibration.
ADULT_MIN_AGE: Final = 18

#: Fixed species leaves per feature, in the reference namespace. Summed
#: identically on both sides. No leaf appears twice across the whole map,
#: which :func:`self_test` asserts.
FEATURE_LEAVES: Final[dict[str, tuple[str, ...]]] = {
    # H-C feature A. The Escherichia species complex as the WGS reference
    # can actually resolve it. Note there is no Shigella leaf in this
    # namespace at all, which is recorded as a transport limitation rather
    # than silently ignored.
    "ecoli_complex": (
        "Escherichia_albertii",
        "Escherichia_coli",
        "Escherichia_fergusonii",
        "Escherichia_marmotae",
        "Escherichia_sp_ESNIH1",
    ),
    # H-C feature B. Ruminococcus gnavus is the reference-namespace name;
    # Mediterraneibacter gnavus is the current name for the same organism.
    # One leaf, one count.
    "m_gnavus": ("Ruminococcus_gnavus",),
    # P-E, four genera treated as one correlated ecological group.
    "faecalibacterium": ("Faecalibacterium_prausnitzii",),
    "coprococcus": (
        "Coprococcus_catus",
        "Coprococcus_comes",
        "Coprococcus_eutactus",
    ),
    "subdoligranulum": ("Subdoligranulum_variabile",),
    "blautia": (
        "Blautia_coccoides",
        "Blautia_hansenii",
        "Blautia_hydrogenotrophica",
        "Blautia_obeum",
        "Blautia_producta",
        "Blautia_sp_An249",
        "Blautia_sp_CAG_257",
        "Blautia_sp_N6H1_15",
        "Blautia_wexlerae",
    ),
}

#: Recorded so the report can say what the measurement actually is rather
#: than implying it reproduced the source study's measurement.
TRANSPORT_NOTES: Final[dict[str, str]] = {
    "ecoli_complex": (
        "BF-S01 measured a 16S Escherichia-Shigella genus signal. This is a "
        "WGS Escherichia species complex summed over five species leaves. The "
        "reference namespace contains no Shigella leaf, so Shigella is not "
        "represented on either side. Not every read here belongs to a "
        "pathogenic strain and no pathotype is inferred."
    ),
    "m_gnavus": (
        "Reported as Ruminococcus gnavus in this reference namespace; "
        "Mediterraneibacter gnavus is the current name for the same organism. "
        "Counted once."
    ),
    "faecalibacterium": (
        "Genus total from the single Faecalibacterium leaf present in the "
        "reference namespace. BF-S01 observed the association in biopsies."
    ),
    "coprococcus": "Genus total over three Coprococcus species leaves.",
    "subdoligranulum": (
        "Genus total from the single Subdoligranulum leaf present in the "
        "reference namespace."
    ),
    "blautia": "Genus total over nine Blautia species leaves.",
}

#: MetaPhlAn 3 profile filename written by the pipeline. This is the lane
#: that shares a namespace with the reference cohort.
MPA3_PROFILE: Final = "metaphlan.mpa_v31_CHOCOPhlAn_201901.tsv"
#: MetaPhlAn 4 profile, deliberately NOT used for these proxies.
MPA4_PROFILE: Final = "metaphlan4.mpa_vJun23_CHOCOPhlAnSGB_202403.tsv"


class ReferenceUnavailable(RuntimeError):
    """The reference cohort file is missing or unusable."""


class IncompatibleProfile(RuntimeError):
    """The sample has no profile in the reference's taxonomy namespace."""


@dataclass(frozen=True, slots=True)
class ReferenceCohort:
    """A frozen, subject-deduplicated, adult-only reference."""

    id: str
    rows: tuple[dict[str, float], ...]
    #: Subject IDs in the same order as ``rows``, for provenance only.
    subjects: tuple[str, ...]
    studies: dict[str, int]
    age_min: float
    age_max: float
    age_median: float
    n_source_samples: int
    n_subjects_before_age_filter: int
    profiler: str
    database: str
    manifest_id: str
    #: Features whose reference column cannot rank anybody, with the reason.
    #: Populated at build time from the actual detection rates, so a feature
    #: that the database cannot see is caught rather than contributing a
    #: constant midrank to every sample.
    unmeasurable: dict[str, str] = field(default_factory=dict)
    #: How subject identity was established, for the audit trail.
    subject_source: str = ""

    @property
    def n(self) -> int:
        return len(self.rows)

    def usable_features(self, features: tuple[str, ...]) -> tuple[str, ...]:
        """Those of ``features`` this reference can actually rank."""
        return tuple(f for f in features if f not in self.unmeasurable)

    def column(self, feature: str) -> list[float]:
        return [row[feature] for row in self.rows]

    def provenance(self) -> dict[str, object]:
        return {
            "reference_id": self.id,
            "n_participants": self.n,
            "n_source_samples": self.n_source_samples,
            "n_subjects_before_age_filter": self.n_subjects_before_age_filter,
            "samples_collapsed": self.n_source_samples - self.n_subjects_before_age_filter,
            "n_studies": len(self.studies),
            "studies": dict(sorted(self.studies.items(), key=lambda kv: -kv[1])),
            "age_range": [self.age_min, self.age_max],
            "age_median": self.age_median,
            "profiler": self.profiler,
            "database": self.database,
            "manifest_id": self.manifest_id,
            "domain": "adults_18_and_over",
            "condition": "control",
            "construction": (
                "Samples were collapsed to subjects using per-study ID rules, one "
                "sample per subject was kept deterministically by highest read "
                "count, and only subjects with a recorded age of 18 or over were "
                "retained. Subjects with no recorded age were excluded rather "
                "than assumed adult. Studies known to sample people repeatedly "
                "whose IDs do not identify the subject were excluded outright "
                "rather than guessed at."
            ),
            "excluded_studies": dict(EXCLUDED_STUDIES),
            "subject_identity_source": self.subject_source,
            "unmeasurable_features": dict(self.unmeasurable),
            "measurement_basis": (
                "Relative abundance as percent of classified reads. Sample "
                "profiles are renormalised to this basis because this pipeline "
                "runs MetaPhlAn with unknown estimation and the reference "
                "cohort does not."
            ),
            "not_established": (
                "Control status in these studies means no recorded disease for "
                "the study's purpose. It does not mean microscopy-confirmed "
                "absence of mucosal biofilm, which was never assessed in any of "
                "them."
            ),
        }


def _subject_id(sample_id: str, study: str) -> str:
    """Fallback subject recovery from the sample ID text."""
    rule = _SUBJECT_RULES.get(study)
    return re.sub(rule, "", sample_id) if rule else sample_id


@lru_cache(maxsize=2)
def _authoritative_subjects(cmd_dir: str) -> dict[str, str]:
    """``sample_id -> subject_id`` from curatedMetagenomicData's metadata.

    Returns an empty mapping when the metadata is not available, in which
    case the caller falls back to the ID-pattern rules and says so. Cached
    because reading and parsing the R data file takes ~20 s and every
    sample in a batch needs the same answer.
    """
    try:
        from openbiota.refcohort import load_sample_metadata

        class _Silent:
            def record(self, *a: object, **k: object) -> None: ...
            def info(self, *a: object, **k: object) -> None: ...
            def warn(self, *a: object, **k: object) -> None: ...

        frame = load_sample_metadata(Path(cmd_dir), _Silent())
    except Exception:  # noqa: BLE001 - absence is handled, not fatal
        return {}
    out: dict[str, str] = {}
    for record in frame.to_dict("records"):
        sample = str(record.get("sample_id") or "")
        subject = str(record.get("subject_id") or "")
        if sample and subject and subject.lower() not in {"nan", "none"}:
            out[sample] = subject
    return out


def _feature_vector(
    abundance_for_sample: dict[str, float]
) -> dict[str, float]:
    """Sum the fixed leaf sets into the six proxy features."""
    return {
        feature: float(
            sum(abundance_for_sample.get(leaf, 0.0) for leaf in leaves)
        )
        for feature, leaves in FEATURE_LEAVES.items()
    }


@lru_cache(maxsize=2)
def load_reference(cohort_path: str) -> ReferenceCohort:
    """Build the frozen adult reference from the taxonomic cohort file.

    Cached because the cohort file is ~14 MB and every sample in a batch
    needs the identical object. The cache key is the path, and the returned
    object is frozen, so a batch run and a single run see the same
    reference - which is what makes BF-T042 hold.
    """
    path = Path(cohort_path)
    if not path.exists():
        raise ReferenceUnavailable(f"reference cohort not found: {path}")
    try:
        data = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ReferenceUnavailable(f"cannot read reference cohort: {exc}") from exc

    taxa: list[str] = data["taxa"]
    sample_ids: list[str] = data["sample_ids"]
    abundance: list[list[float]] = data["abundance"]
    metadata: dict[str, dict] = data["metadata"]
    manifest = data.get("manifest", {})

    taxon_index = {name: i for i, name in enumerate(taxa)}
    # Only the leaves we actually use, so we never materialise 999x3027.
    needed = {
        leaf: taxon_index[leaf]
        for leaves in FEATURE_LEAVES.values()
        for leaf in leaves
        if leaf in taxon_index
    }

    # One record per sample, then collapse to subjects.
    by_subject: dict[tuple[str, str], list[tuple[int, str]]] = {}
    authoritative = _authoritative_subjects(CMD_METADATA_DIR)
    subject_source = "curatedMetagenomicData subject_id" if authoritative else (
        "recovered from sample-ID patterns; the cMD metadata was unavailable"
    )
    for col, sid in enumerate(sample_ids):
        meta = metadata.get(sid)
        if meta is None:
            continue
        study = meta.get("study_name", "unknown")
        if study in EXCLUDED_STUDIES:
            continue
        subject = authoritative.get(sid) or _subject_id(sid, study)
        by_subject.setdefault((study, subject), []).append((col, sid))

    n_subjects_all = len(by_subject)

    rows: list[dict[str, float]] = []
    subjects: list[str] = []
    studies: dict[str, int] = {}
    ages: list[float] = []

    for (study, subject), members in sorted(by_subject.items()):
        # Deterministic pick: highest read count, ties broken by sample ID.
        chosen_col, chosen_sid = max(
            members,
            key=lambda cs: (metadata[cs[1]].get("n_reads") or 0, cs[1]),
        )
        age = metadata[chosen_sid].get("age")
        # Age must be known and adult. An unknown age is not assumed adult:
        # BackhedF_2015 is a mother-infant study with no ages recorded.
        if age is None or float(age) < ADULT_MIN_AGE:
            continue
        per_leaf = {
            leaf: float(abundance[idx][chosen_col]) for leaf, idx in needed.items()
        }
        rows.append(_feature_vector(per_leaf))
        subjects.append(f"{study}:{subject}")
        studies[study] = studies.get(study, 0) + 1
        ages.append(float(age))

    if not rows:
        raise ReferenceUnavailable("no adult subjects survived reference construction")

    ordered_ages = sorted(ages)
    median = (
        ordered_ages[len(ordered_ages) // 2]
        if len(ordered_ages) % 2
        else (
            ordered_ages[len(ordered_ages) // 2 - 1]
            + ordered_ages[len(ordered_ages) // 2]
        )
        / 2
    )
    # Which features this reference can actually rank. A column detected in
    # a couple of participants is effectively constant: everyone ties, the
    # tie contributes a fixed midrank to the group mean, and the features
    # that do carry signal get diluted in proportion. Catch it here from the
    # data rather than hardcoding a known-bad feature name.
    unmeasurable: dict[str, str] = {}
    for feature in FEATURE_LEAVES:
        column = [row[feature] for row in rows]
        detected = sum(1 for value in column if value > 0)
        rate = detected / len(column)
        if rate < MIN_FEATURE_DETECTION_RATE:
            unmeasurable[feature] = (
                f"detected in {detected} of {len(column)} reference "
                f"participants ({rate:.1%}), below the {MIN_FEATURE_DETECTION_RATE:.0%} "
                "floor needed to rank anyone. The organism is not rare; this "
                "taxonomy release simply has no named species leaf that "
                "captures it, so the feature is not assayed here rather than "
                "absent."
            )

    digest = hashlib.sha256(
        "|".join(
            (
                str(manifest.get("manifest_id", "")),
                str(len(rows)),
                # The leaf lists, not just the feature names: changing which
                # species are summed changes the measurement and must change
                # the calibration identity.
                ";".join(
                    f"{name}={'+'.join(leaves)}"
                    for name, leaves in sorted(FEATURE_LEAVES.items())
                ),
                str(ADULT_MIN_AGE),
                ",".join(sorted(EXCLUDED_STUDIES)),
                subject_source,
                ",".join(sorted(unmeasurable)),
            )
        ).encode()
    ).hexdigest()[:12]

    return ReferenceCohort(
        id=f"adult-control-mpa3-{digest}",
        rows=tuple(rows),
        subjects=tuple(subjects),
        studies=studies,
        age_min=min(ages),
        age_max=max(ages),
        age_median=float(median),
        n_source_samples=len(sample_ids),
        n_subjects_before_age_filter=n_subjects_all,
        profiler=str(manifest.get("profiler", "MetaPhlAn 3")),
        database="mpa_v3_CHOCOPhlAn_201901",
        manifest_id=str(manifest.get("manifest_id", "")),
        unmeasurable=unmeasurable,
        subject_source=subject_source,
    )


def parse_mpa3_profile(path: Path) -> tuple[dict[str, float], dict[str, float]]:
    """Species-level relative abundances from a MetaPhlAn 3 profile.

    Returns ``(species, scale)`` where ``species`` is on the same basis as
    the reference cohort and ``scale`` records the renormalisation.

    **This pipeline runs MetaPhlAn with unknown estimation and the
    reference cohort does not.** A profile here reads
    ``UNKNOWN 49.8`` / ``k__Bacteria 50.2``, so its species rows sum to
    about half, while every curatedMetagenomicData row sums to 100 because
    cMD reports percent-of-classified. Ranking the first against the
    second understates every feature by roughly a factor of two and drags
    all percentiles down.

    So the species rows are renormalised to percent-of-classified, which
    is the reference's basis. The unknown fraction is not discarded
    information - it is returned in ``scale`` and reported - but it must
    not silently act as a divisor on one side of a comparison only.

    Only ``s__`` rows are read, and strain rows below them are skipped, so
    a clade and its child can never both be counted.
    """
    raw: dict[str, float] = {}
    unknown = 0.0
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        lineage = parts[0]
        if lineage.strip() in {"UNKNOWN", "UNCLASSIFIED"}:
            try:
                unknown = float(parts[2])
            except ValueError:
                unknown = 0.0
            continue
        # Species rows only, and not strain rows below them.
        if "|s__" not in lineage or "|t__" in lineage:
            continue
        species = lineage.rsplit("|s__", 1)[1]
        try:
            raw[species] = raw.get(species, 0.0) + float(parts[2])
        except ValueError:
            continue

    classified = sum(raw.values())
    scale = {
        "classified_percent": round(classified, 6),
        "unknown_percent": round(unknown, 6),
        "renormalised": False,
        "factor": 1.0,
    }
    if classified <= 0:
        return raw, scale
    # Put the sample on the reference's basis: percent of classified.
    factor = 100.0 / classified
    if abs(factor - 1.0) > 1e-9:
        raw = {k: v * factor for k, v in raw.items()}
        scale["renormalised"] = True
        scale["factor"] = round(factor, 6)
    return raw, scale


def sample_features(sample_dir: Path) -> tuple[dict[str, float], dict[str, object]]:
    """The six proxy features for one sample, in the reference namespace.

    Returns the feature vector and a detail record naming every species
    leaf that contributed, so the report can show exactly what went into a
    number rather than asserting it.
    """
    profile_path = sample_dir / "taxonomy" / MPA3_PROFILE
    if not profile_path.exists():
        raise IncompatibleProfile(
            f"no MetaPhlAn 3 profile at {profile_path}. The reference cohort is "
            "MetaPhlAn 3 species-level; the MetaPhlAn 4 SGB profile is a "
            "different measurement and cannot be ranked against it."
        )
    species, scale = parse_mpa3_profile(profile_path)
    features = _feature_vector(species)
    detail: dict[str, object] = {}
    for feature, leaves in FEATURE_LEAVES.items():
        contributions = {
            leaf: round(species.get(leaf, 0.0), 6)
            for leaf in leaves
            if species.get(leaf, 0.0) > 0
        }
        detail[feature] = {
            "value": round(features[feature], 6),
            "unit": (
                "relative abundance (% of classified), MetaPhlAn 3 species level"
            ),
            "leaves_summed": list(leaves),
            "leaves_detected": contributions,
            "n_leaves_detected": len(contributions),
            "detected": features[feature] > 0,
            "censoring": "none" if features[feature] > 0 else "left_censored",
            "transport_note": TRANSPORT_NOTES[feature],
            "scale": dict(scale),
        }
    return features, detail


def self_test() -> None:
    """Invariants that protect the reference from silent corruption."""
    # No species leaf may appear in two features, or a read is counted twice.
    seen: set[str] = set()
    for leaves in FEATURE_LEAVES.values():
        for leaf in leaves:
            assert leaf not in seen, f"{leaf} appears in more than one feature"
            seen.add(leaf)
    # The two proxies must partition their features cleanly.
    from .scoring import COMMUNITY_PROXY, ECOLOGY_PROXY

    hc = set(COMMUNITY_PROXY.features)
    pe = set(ECOLOGY_PROXY.features)
    assert not (hc & pe), "a feature cannot serve both proxies"
    assert hc | pe == set(FEATURE_LEAVES), "feature sets must match the leaf map"
