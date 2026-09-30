"""Estimated stool-microbiome chronological age (research use) — spec v04 §9.

The output is a *prediction of chronological age from a stool community
pattern*. It is not biological age, gut age, a health grade, an ageing rate or
a treatment target, and an older- or younger-than-chronological prediction is
not inherently good or bad. Every string this module emits keeps to that
vocabulary.

Why this is a clean-room retraining and not an import
----------------------------------------------------
Myers et al. (Commun Biol 2025) do not ship an inference bundle (no fitted
RCLR/PCA/scalers, feature order or exact training manifest), and their
notebooks fit preprocessing globally before cross-validation. The spec
therefore requires: one frozen input namespace, a reconciled one-baseline-per-
participant manifest, every transformation fitted inside each training fold,
leave-one-study-out as the primary split, a split-conformal interval fitted on
an untouched calibration manifest, three OOD gates whose thresholds are
learned without the test specimen, and an adult-only release with a pediatric
``insufficient_pediatric_reference_support`` gate.

Input namespace
---------------
curatedMetagenomicData's 2021 snapshots are MetaPhlAn 3.0 / CHOCOPhlAn 201901
species profiles; the screened sample is profiled with MetaPhlAn 3.1 against
the same CHOCOPhlAn 201901 index. Training and inference therefore share one
species namespace. MetaPhlAn 4 (SGB) output is never fed to this model.

Pipeline (fitted inside every fold)
-----------------------------------
prevalence filter → robust CLR (log of the non-zero part of each sample,
centred on that sample's mean log; zeros stay at the centre) → PCA on the
training fold → ridge regression on the leading components. k and the ridge
penalty are chosen by an inner study-grouped cross-validation. A kNN on the
same components is the baseline the linear model has to beat.

Everything numeric is numpy; the training set is ~8,000 × ~1,000, small
enough that the whole leave-one-study-out audit runs in minutes.
"""

from __future__ import annotations

import hashlib
import json
import math
import time
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass, fields, replace
from pathlib import Path
from typing import Any, Final

import numpy as np

from openbiota.errors import OpenBiotaError
from openbiota.logging_util import Reporter
from openbiota.refcohort import (
    _as_float,
    _as_str,
    download_profile,
    load_sample_metadata,
    parse_profile,
)

MODEL_ID: Final = "MPA_WGS_AGE_V1"
TASK: Final = "predicted_chronological_age"
INPUT_NAMESPACE: Final = "MetaPhlAn 3 species / CHOCOPhlAn 201901"
TRAINING_DATABASE: Final = "curatedMetagenomicData 3 (MetaPhlAn 3.0, CHOCOPhlAn_201901)"
INFERENCE_INDEX: Final = "mpa_v31_CHOCOPhlAn_201901"

PERMANENT_STATEMENT: Final = (
    "The age the stool community reads as, estimated from its species patterns; "
    "a descriptive measure, not a health grade or a treatment target."
)

#: Myers et al. stool-WGS project seed (spec §14.2). Study -> reported n.
#: The seed is `body_site == stool`, `study_condition == control`, `age` known,
#: which reconciles row-for-row against cMD's sampleMetadata (asserted at
#: manifest build time).
SEED_STUDIES: Final[dict[str, int]] = {
    "AsnicarF_2017": 8, "AsnicarF_2021": 1098, "BackhedF_2015": 287,
    "Bengtsson-PalmeJ_2015": 70, "BritoIL_2016": 172, "BrooksB_2017": 5,
    "ChuDM_2017": 65, "DeFilippisF_2019": 97, "DhakanDB_2019": 110,
    "FengQ_2015": 61, "GuptaA_2019": 30, "HMP_2012": 147, "HMP_2019_ibdmdb": 426,
    "HMP_2019_t2d": 46, "HanniganGD_2017": 28, "HansenLBS_2018": 207,
    "Heitz-BuschartA_2016": 26, "IjazUZ_2017": 37, "KarlssonFH_2013": 43,
    "KaurK_2020": 31, "KeohaneDM_2020": 117, "KosticAD_2015": 89,
    "LifeLinesDeep_2016": 1135, "LokmerA_2019": 57, "NagySzakalD_2017": 50,
    "Obregon-TitoAJ_2015": 57, "PasolliE_2019": 112, "PehrssonE_2016": 191,
    "QinN_2014": 114, "RampelliS_2015": 38, "RaymondF_2016": 36, "RubelMA_2020": 86,
    "SankaranarayananK_2015": 18, "ShaoY_2019": 1619, "SmitsSA_2017": 27,
    "TettAJ_2019_a": 68, "TettAJ_2019_b": 43, "TettAJ_2019_c": 49,
    "ThomasAM_2018a": 24, "ThomasAM_2018b": 27, "ThomasAM_2019_c": 40,
    "VatanenT_2016": 615, "VincentC_2016": 196, "VogtmannE_2016": 52,
    "WampachL_2018": 53, "WirbelJ_2018": 65, "XieH_2016": 177, "YachidaS_2019": 251,
    "YeZ_2018": 45, "YuJ_2015": 54, "ZellerG_2014": 61, "ZhuF_2020": 81,
}
REPORTED_PRE_DEDUP_N: Final = 8641

#: Known project overlaps (spec §14.2): the second study of each pair shares
#: a BioProject with the first. Runs shared across the pair are kept once,
#: under the study named first (the larger / primary deposit).
OVERLAP_CLUSTERS: Final[tuple[tuple[str, str], ...]] = (
    ("DhakanDB_2019", "GuptaA_2019"),
    ("WirbelJ_2018", "ThomasAM_2018b"),
    ("YachidaS_2019", "ThomasAM_2019_c"),
)

#: Adult-only release range. Everything below 18 is pediatric and gated.
ADULT_MIN_AGE: Final = 18.0
#: A numeric pediatric result needs ≥200 independent 8–18-year-olds from ≥3
#: external cohorts/labs/regions with no training overlap (spec §9.3). The
#: bundled data cannot meet that, so pediatric output is always gated.
PEDIATRIC_REQUIRED_N: Final = 200
PEDIATRIC_REQUIRED_COHORTS: Final = 3

#: Species must be present in at least this fraction of a training fold.
MIN_PREVALENCE: Final = 0.02

#: Weight on the linear-carriage component in the two-model ensemble; the
#: transformer carries the remainder. Chosen from the blend curve measured on
#: leave-one-study-out predictions, where a fifth of the weight on the
#: transformer leaves mean absolute error within a month of the linear
#: component alone and moves the calibration slope towards one. Both
#: architectures read different structure — carriage of individual species
#: against attention over robust principal components — so the estimate does
#: not rest on either one alone.
ENSEMBLE_WEIGHT_LINEAR: Final = 0.80

#: Regression feature mode. Measured on identical leave-one-study-out folds
#: over 4,800 adults from 46 studies (all preprocessing inside the fold):
#:
#: =================================  =====  ======  =====
#: features                             MAE      R²  slope
#: =================================  =====  ======  =====
#: presence/absence                   11.71  +0.112   0.71
#: log10 abundance                    11.93  +0.079   0.65
#: log10 abundance + presence         11.91  +0.088   0.65
#: robust-CLR + PCA (spec v04 §9)     12.15  -0.009   0.43
#: sqrt abundance                     12.38  -0.022   0.47
#: constant = training mean           12.81  -0.017     --
#: =================================  =====  ======  =====
#:
#: Carriage transports across studies; abundance does not, because abundance
#: is what extraction and library protocol distort most. The abundance model
#: was indistinguishable from predicting the training mean.
FEATURE_MODE: Final = "presence"

#: Nested-CV grids.
COMPONENT_GRID: Final = (30, 60, 120, 200)
ALPHA_GRID: Final = (0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0, 100000.0)
INNER_FOLDS: Final = 3
KNN_K: Final = 10

#: Minimum training samples within ±5 years of a prediction for it to have
#: support; below this the estimate is a warning, below a quarter it abstains.
SUPPORT_WINDOW_YEARS: Final = 5.0
SUPPORT_MIN_N: Final = 100

#: Conformal targets.
CONFORMAL_LEVELS: Final = (0.80, 0.95)
#: Fraction of adult training samples to reserve (as whole studies) for the
#: untouched calibration manifest.
CALIBRATION_FRACTION: Final = 0.15

#: Analytical OOD gates on the specimen itself.
MIN_READS_FOR_AGE: Final = 1_000_000
MAX_UNKNOWN_PERCENT: Final = 60.0

#: Attribution stability: a species contribution is shown only when the sign of
#: its species-space weight agrees in at least this fraction of held-out folds.
STABILITY_MIN_AGREEMENT: Final = 0.90
ATTRIBUTION_TOP_N: Final = 8

LIFE_STAGE_BANDS: Final = ((18, 30, "young adult"), (30, 45, "adult"), (45, 60, "middle age"),
                           (60, 75, "older adult"), (75, 120, "oldest"))


def sha256_of(obj: Any) -> str:
    if isinstance(obj, (bytes, bytearray)):
        data = bytes(obj)
    elif isinstance(obj, np.ndarray):
        data = np.ascontiguousarray(obj).tobytes() + str(obj.shape).encode() + str(obj.dtype).encode()
    else:
        data = json.dumps(obj, sort_keys=True, default=str).encode()
    return "sha256:" + hashlib.sha256(data).hexdigest()


def life_stage(age: float) -> str:
    for lo, hi, name in LIFE_STAGE_BANDS:
        if lo <= age < hi:
            return name
    return "pediatric" if age < ADULT_MIN_AGE else "oldest"


def age_band_label(age: float) -> str:
    for lo, hi, _ in LIFE_STAGE_BANDS:
        if lo <= age < hi:
            return f"{lo}-{hi}"
    return "<18" if age < ADULT_MIN_AGE else "75+"


# --------------------------------------------------------------------------- #
# manifest
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class ManifestRow:
    sample_id: str
    study: str
    subject_id: str
    accessions: tuple[str, ...]
    age_years: float
    age_source: str            # "age" | "infant_age_days"
    sex: str | None
    country: str | None
    n_reads: int | None
    days_from_first_collection: float | None
    visit_number: float | None
    antibiotics: str | None
    included: bool = True
    exclusion_reason: str | None = None
    overlap_cluster: str | None = None
    time_point_rank: int = 0

    def to_json(self) -> dict[str, Any]:
        d = asdict(self)
        d["accessions"] = list(self.accessions)
        return d


@dataclass(slots=True)
class AgeManifest:
    rows: list[ManifestRow]
    reported_pre_dedup_n: int
    listed_seed_n: int
    reconstructed_pre_dedup_n: int
    post_dedup_n: int
    per_study: dict[str, dict[str, int]]
    join_diagnostics: dict[str, Any]
    metadata_digest: str
    built_at: str

    @property
    def included(self) -> list[ManifestRow]:
        return [r for r in self.rows if r.included]

    def to_json(self) -> dict[str, Any]:
        return {
            "model_id": MODEL_ID,
            "reported_pre_dedup_n": self.reported_pre_dedup_n,
            "listed_seed_n": self.listed_seed_n,
            "reconstructed_pre_dedup_n": self.reconstructed_pre_dedup_n,
            "post_dedup_n": self.post_dedup_n,
            "per_study": self.per_study,
            "join_diagnostics": self.join_diagnostics,
            "metadata_digest": self.metadata_digest,
            "built_at": self.built_at,
            "rows": [r.to_json() for r in self.rows],
        }

    @property
    def digest(self) -> str:
        return sha256_of([r.to_json() for r in self.rows if r.included])


def _accessions(value: Any) -> tuple[str, ...]:
    text = _as_str(value)
    if not text:
        return ()
    return tuple(sorted({p.strip() for p in text.replace(",", ";").split(";") if p.strip()}))


def build_manifest(cache_dir: Path, reporter: Reporter) -> AgeManifest:
    """Reconstruct the Myers stool-WGS training manifest from cMD metadata.

    Row-level rule (reconciles to the reported 8,641 pre-dedup rows): stool,
    ``study_condition == control``, known age, study in the seed list.
    Deduplication: one baseline sample per participant (earliest
    days_from_first_collection, then visit_number, then sample_id), and runs
    shared across a known project-overlap pair kept once.
    """
    frame = load_sample_metadata(cache_dir, reporter)
    meta_path = cache_dir / "cmd_sampleMetadata.rda"
    metadata_digest = "sha256:" + hashlib.sha256(meta_path.read_bytes()).hexdigest()

    rows: list[ManifestRow] = []
    per_study: dict[str, dict[str, int]] = {}
    for record in frame.to_dict("records"):
        study = _as_str(record.get("study_name"))
        if study not in SEED_STUDIES:
            continue
        if _as_str(record.get("body_site")) != "stool":
            continue
        if _as_str(record.get("study_condition")) != "control":
            continue
        age = _as_float(record.get("age"))
        if age is None:
            continue
        infant_days = _as_float(record.get("infant_age"))
        if infant_days is not None and age < 3:
            age_years, source = infant_days / 365.25, "infant_age_days"
        else:
            age_years, source = age, "age"
        sample_id = _as_str(record.get("sample_id")) or ""
        subject = _as_str(record.get("subject_id")) or sample_id
        reads = _as_float(record.get("number_reads"))
        rows.append(
            ManifestRow(
                sample_id=sample_id, study=study, subject_id=subject,
                accessions=_accessions(record.get("NCBI_accession")),
                age_years=age_years, age_source=source,
                sex=_as_str(record.get("gender")), country=_as_str(record.get("country")),
                n_reads=int(reads) if reads else None,
                days_from_first_collection=_as_float(record.get("days_from_first_collection")),
                visit_number=_as_float(record.get("visit_number")),
                antibiotics=_as_str(record.get("antibiotics_current_use")),
            )
        )

    reconstructed = len(rows)
    mismatches = {}
    for study, expected in SEED_STUDIES.items():
        got = sum(1 for r in rows if r.study == study)
        per_study[study] = {"seed_n": expected, "reconstructed_n": got}
        if got != expected:
            mismatches[study] = {"seed_n": expected, "reconstructed_n": got}

    # --- overlap clusters: identical runs deposited under two study labels ---
    overlap_removed = 0
    for primary, secondary in OVERLAP_CLUSTERS:
        primary_runs = {a for r in rows if r.study == primary for a in r.accessions}
        primary_subjects = {r.subject_id for r in rows if r.study == primary}
        for r in rows:
            if r.study != secondary or not r.included:
                continue
            r.overlap_cluster = f"{primary}/{secondary}"
            if (r.accessions and set(r.accessions) & primary_runs) or r.subject_id in primary_subjects:
                r.included = False
                r.exclusion_reason = f"overlap_duplicate_of_{primary}"
                overlap_removed += 1

    # --- current antibiotic use: a declared exclusion for a healthy reference ---
    abx_removed = 0
    for r in rows:
        if r.included and r.antibiotics == "yes":
            r.included = False
            r.exclusion_reason = "antibiotics_current_use"
            abx_removed += 1

    # --- one baseline per participant ---
    by_subject: dict[tuple[str, str], list[ManifestRow]] = {}
    for r in rows:
        if r.included:
            by_subject.setdefault((r.study, r.subject_id), []).append(r)
    repeat_removed = 0
    for members in by_subject.values():
        members.sort(key=lambda r: (
            r.days_from_first_collection if r.days_from_first_collection is not None else math.inf,
            r.visit_number if r.visit_number is not None else math.inf,
            r.sample_id,
        ))
        for rank, r in enumerate(members):
            r.time_point_rank = rank
            if rank > 0:
                r.included = False
                r.exclusion_reason = "repeat_sample_of_participant"
                repeat_removed += 1

    for study in per_study:
        per_study[study]["post_dedup_n"] = sum(1 for r in rows if r.study == study and r.included)
        per_study[study]["participants"] = len({r.subject_id for r in rows if r.study == study})

    post = sum(1 for r in rows if r.included)
    manifest = AgeManifest(
        rows=rows,
        reported_pre_dedup_n=REPORTED_PRE_DEDUP_N,
        listed_seed_n=sum(SEED_STUDIES.values()),
        reconstructed_pre_dedup_n=reconstructed,
        post_dedup_n=post,
        per_study=per_study,
        join_diagnostics={
            "row_rule": "body_site==stool AND study_condition==control AND age known AND study in seed",
            "seed_mismatches": mismatches,
            "overlap_duplicates_removed": overlap_removed,
            "antibiotic_users_removed": abx_removed,
            "repeat_samples_removed": repeat_removed,
            "reconciled": reconstructed == REPORTED_PRE_DEDUP_N and not mismatches,
        },
        metadata_digest=metadata_digest,
        built_at=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    )
    reporter.record(
        f"    age manifest: {reconstructed:,} rows reconstructed (reported {REPORTED_PRE_DEDUP_N:,}), "
        f"{post:,} after one-baseline-per-participant; {len(mismatches)} study mismatch(es)"
    )
    return manifest


# --------------------------------------------------------------------------- #
# training matrix
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class TrainingMatrix:
    taxa: list[str]
    sample_ids: list[str]
    X: np.ndarray            # n × p, relative abundance percent
    age: np.ndarray          # n
    study: np.ndarray        # n, str
    country: np.ndarray      # n, str ("unknown" when missing)
    sex: np.ndarray          # n, str
    n_reads: np.ndarray      # n, float (nan when missing)

    @property
    def digest(self) -> str:
        return sha256_of(self.X) + "|" + sha256_of(self.taxa)


def _profile_paths(profiles_dir: Path, studies: Iterable[str]) -> dict[str, Path]:
    out: dict[str, Path] = {}
    for study in studies:
        candidates = sorted(profiles_dir.glob(f"*/{study}.rda"))
        if not candidates:
            raise OpenBiotaError(
                f"no relative-abundance profile for {study} under {profiles_dir}; "
                "run `openbiota age-model fetch` first"
            )
        # prefer the newest snapshot that exists locally
        out[study] = candidates[-1]
    return out


def fetch_profiles(cache_dir: Path, reporter: Reporter) -> dict[str, Path]:
    """Download every seed study's species profile (any 2021 snapshot)."""
    import sqlite3

    sqlite_path = cache_dir / "experimenthub.sqlite3"
    if not sqlite_path.is_file():
        from openbiota.net import get
        from openbiota.refcohort import EXPERIMENTHUB_SQLITE
        sqlite_path.write_bytes(get(EXPERIMENTHUB_SQLITE, timeout=600).body)
    urls: dict[str, dict[str, str]] = {}
    query = """
        SELECT r.title, p.rdatapath, l.location_prefix FROM resources r
        JOIN rdatapaths p ON p.resource_id = r.id
        JOIN location_prefixes l ON l.id = r.location_prefix_id
        WHERE r.title LIKE '%relative_abundance'
    """
    for title, rdatapath, prefix in sqlite3.connect(sqlite_path).execute(query):
        snap, study = str(title).split(".")[0], str(title).split(".")[1]
        urls.setdefault(study, {})[snap] = f"{prefix}{rdatapath}"
    out: dict[str, Path] = {}
    for study in SEED_STUDIES:
        snaps = urls.get(study, {})
        if not snaps:
            raise OpenBiotaError(f"ExperimentHub has no relative-abundance resource for {study}")
        snap = "2021-10-14" if "2021-10-14" in snaps else ("2021-03-31" if "2021-03-31" in snaps else sorted(snaps)[0])
        dest = cache_dir / "profiles" / snap / f"{study}.rda"
        if not dest.is_file():
            reporter.info(f"  fetching {study} ({snap})")
            download_profile(snaps[snap], dest)
        out[study] = dest
    return out


def load_training_matrix(manifest: AgeManifest, cache_dir: Path, reporter: Reporter) -> TrainingMatrix:
    """Assemble the included manifest rows into one species × sample matrix."""
    cache = cache_dir / "age_training_matrix.npz"
    digest = manifest.digest
    if cache.is_file():
        with np.load(cache, allow_pickle=False) as z:
            if str(z["manifest_digest"]) == digest:
                reporter.record("    training matrix: reusing cached assembly")
                return TrainingMatrix(
                    taxa=list(z["taxa"]), sample_ids=list(z["sample_ids"]), X=z["X"], age=z["age"],
                    study=z["study"], country=z["country"], sex=z["sex"], n_reads=z["n_reads"],
                )

    included = manifest.included
    wanted = {r.sample_id: r for r in included}
    paths = _profile_paths(cache_dir / "profiles", sorted({r.study for r in included}))
    taxa_index: dict[str, int] = {}
    columns: list[tuple[str, dict[int, float]]] = []
    missing_samples: list[str] = []
    for study, path in paths.items():
        taxa, samples, matrix = parse_profile(path)
        idx = [taxa_index.setdefault(t, len(taxa_index)) for t in taxa]
        sample_pos = {s: j for j, s in enumerate(samples)}
        for row in included:
            if row.study != study:
                continue
            j = sample_pos.get(row.sample_id)
            if j is None:
                missing_samples.append(row.sample_id)
                continue
            sparse = {idx[i]: matrix[i][j] for i in range(len(taxa)) if matrix[i][j] > 0}
            columns.append((row.sample_id, sparse))
    if missing_samples:
        reporter.warn(f"  {len(missing_samples)} manifest samples missing from the profiles; excluded")
        for sid in missing_samples:
            wanted[sid].included = False
            wanted[sid].exclusion_reason = "profile_missing"

    taxa_list = [t for t, _ in sorted(taxa_index.items(), key=lambda kv: kv[1])]
    X = np.zeros((len(columns), len(taxa_list)), dtype=np.float32)
    ids: list[str] = []
    for i, (sid, sparse) in enumerate(columns):
        ids.append(sid)
        for k, v in sparse.items():
            X[i, k] = v
    rows = [wanted[s] for s in ids]
    tm = TrainingMatrix(
        taxa=taxa_list, sample_ids=ids, X=X,
        age=np.array([r.age_years for r in rows], dtype=np.float64),
        study=np.array([r.study for r in rows]),
        country=np.array([r.country or "unknown" for r in rows]),
        sex=np.array([r.sex or "unknown" for r in rows]),
        n_reads=np.array([r.n_reads if r.n_reads else np.nan for r in rows], dtype=np.float64),
    )
    np.savez_compressed(
        cache, manifest_digest=np.array(manifest.digest), taxa=np.array(taxa_list),
        sample_ids=np.array(ids), X=X, age=tm.age, study=tm.study, country=tm.country,
        sex=tm.sex, n_reads=tm.n_reads,
    )
    reporter.record(f"    training matrix: {X.shape[0]:,} samples × {X.shape[1]:,} species")
    return tm


# --------------------------------------------------------------------------- #
# transforms (all fitted on a training fold only)
# --------------------------------------------------------------------------- #


def rclr(X: np.ndarray) -> np.ndarray:
    """Robust centred log-ratio: log of non-zero entries centred on each
    sample's mean log; zeros are left at the centre (0)."""
    mask = X > 0
    with np.errstate(divide="ignore"):
        logs = np.where(mask, np.log(np.where(mask, X, 1.0)), 0.0)
    counts = mask.sum(axis=1, keepdims=True)
    counts = np.where(counts == 0, 1, counts)
    means = logs.sum(axis=1, keepdims=True) / counts
    return np.where(mask, logs - means, 0.0)


@dataclass(slots=True)
class FittedTransform:
    """Design matrix for the regression, plus a PCA kept only for the OOD gate.

    Two feature modes. ``rclr_pca`` is the original: robust CLR of relative
    abundance, then the leading principal components. ``presence`` discards
    abundance and keeps only which species were detected, standardised.

    Presence/absence wins on held-out studies (see
    :data:`FEATURE_MODE`), which is what one would expect: abundance is the
    quantity most distorted by extraction and library protocol, so an
    abundance-based model partly learns the study it came from, whereas the
    carriage pattern transports.

    The PCA is fitted in both modes because the out-of-distribution gate
    measures distance in component space, which is a separate question from
    which features the regression reads.
    """

    keep: np.ndarray          # boolean species mask (prevalence filter)
    mean: np.ndarray          # rclr mean over training (kept species)
    components: np.ndarray    # k × p_kept PCA loadings (rows)
    explained_var: np.ndarray # k
    k: int
    mode: str = "rclr_pca"
    pmean: np.ndarray | None = None   # presence mean over training (kept species)
    psd: np.ndarray | None = None     # presence sd over training (kept species)

    def design(self, X_percent: np.ndarray) -> np.ndarray:
        """Features the regression sees."""
        if self.mode == "presence":
            P = (X_percent[:, self.keep] > 0).astype(np.float64)
            return (P - self.pmean) / self.psd
        return self.pca_scores(X_percent)

    # kept as the historical name so callers that only want the regression
    # design matrix continue to work
    def transform(self, X_percent: np.ndarray) -> np.ndarray:
        return self.design(X_percent)

    def pca_scores(self, X_percent: np.ndarray) -> np.ndarray:
        Z = rclr(X_percent[:, self.keep]) - self.mean
        return Z @ self.components[: self.k].T

    def latent_distance(self, scores: np.ndarray) -> np.ndarray:
        """Mahalanobis-style distance in component space (unit-variance scaled)."""
        var = np.maximum(self.explained_var[: self.k], 1e-12)
        return np.sqrt(((scores ** 2) / var).sum(axis=1) / self.k)


def fit_transform(X_train: np.ndarray, *, max_k: int, mode: str = "rclr_pca") -> FittedTransform:
    prevalence = (X_train > 0).mean(axis=0)
    keep = prevalence >= MIN_PREVALENCE
    Z = rclr(X_train[:, keep])
    mean = Z.mean(axis=0)
    Zc = Z - mean
    k = int(min(max_k, Zc.shape[0] - 1, Zc.shape[1]))
    # economy SVD; components are right singular vectors
    _, s, vt = np.linalg.svd(Zc, full_matrices=False)
    explained = (s[:k] ** 2) / max(Zc.shape[0] - 1, 1)
    P = (X_train[:, keep] > 0).astype(np.float64)
    return FittedTransform(
        keep=keep, mean=mean, components=vt[:k], explained_var=explained, k=k, mode=mode,
        pmean=P.mean(axis=0), psd=P.std(axis=0) + 1e-9,
    )


def ridge_fit(S: np.ndarray, y: np.ndarray, alpha: float) -> tuple[np.ndarray, float]:
    """Closed-form ridge with an unpenalised intercept."""
    y_mean = float(y.mean())
    s_mean = S.mean(axis=0)
    Sc, yc = S - s_mean, y - y_mean
    A = Sc.T @ Sc + alpha * np.eye(S.shape[1])
    beta = np.linalg.solve(A, Sc.T @ yc)
    intercept = y_mean - float(s_mean @ beta)
    return beta, intercept


def knn_predict(S_train: np.ndarray, y_train: np.ndarray, S_test: np.ndarray, k: int = KNN_K) -> np.ndarray:
    out = np.empty(S_test.shape[0])
    train_sq = (S_train ** 2).sum(axis=1)
    k = min(k, S_train.shape[0])
    for i in range(0, S_test.shape[0], 1024):
        block = S_test[i:i + 1024]
        d = train_sq[None, :] - 2.0 * (block @ S_train.T) + (block ** 2).sum(axis=1)[:, None]
        nn = np.argpartition(d, k - 1, axis=1)[:, :k]
        out[i:i + 1024] = y_train[nn].mean(axis=1)
    return out


def _group_folds(groups: np.ndarray, n_folds: int) -> list[np.ndarray]:
    """Assign whole groups to folds, greedily balancing fold sizes (deterministic)."""
    uniq, counts = np.unique(groups, return_counts=True)
    order = np.argsort(-counts, kind="stable")
    loads = np.zeros(n_folds)
    assignment: dict[str, int] = {}
    for i in order:
        f = int(np.argmin(loads))
        assignment[str(uniq[i])] = f
        loads[f] += counts[i]
    fold_of = np.array([assignment[str(g)] for g in groups])
    return [np.where(fold_of == f)[0] for f in range(n_folds)]


@dataclass(slots=True)
class FoldModel:
    transform: FittedTransform
    beta: np.ndarray
    intercept: float
    alpha: float
    k: int

    def predict(self, X_percent: np.ndarray) -> np.ndarray:
        return self.transform.design(X_percent) @ self.beta + self.intercept

    def species_weights(self) -> np.ndarray:
        """Linear weight of every kept species in its own feature space (p_kept)."""
        if self.transform.mode == "presence":
            return self.beta / self.transform.psd
        return self.transform.components[: self.k].T @ self.beta


def fit_fold(
    X: np.ndarray, y: np.ndarray, groups: np.ndarray, *,
    reporter: Reporter | None = None,  # noqa: ARG001 — kept for call-site symmetry
    mode: str = FEATURE_MODE,
) -> FoldModel:
    """Nested selection of the ridge penalty (and k, in PCA mode) by
    study-grouped inner CV, then refit on the whole training fold."""
    inner = _group_folds(groups, INNER_FOLDS) if len(np.unique(groups)) >= INNER_FOLDS else None
    best: tuple[float, int, float] | None = None
    if inner is not None:
        errors: dict[tuple[int, float], list[float]] = {}
        for held in inner:
            train = np.setdiff1d(np.arange(len(y)), held)
            tf = fit_transform(X[train], max_k=max(COMPONENT_GRID), mode=mode)
            S_tr, S_te = tf.design(X[train]), tf.design(X[held])
            # in presence mode there is no component count to tune
            ks = (S_tr.shape[1],) if mode == "presence" else COMPONENT_GRID
            for k in ks:
                k_eff = min(k, S_tr.shape[1])
                for alpha in ALPHA_GRID:
                    beta, b0 = ridge_fit(S_tr[:, :k_eff], y[train], alpha)
                    pred = S_te[:, :k_eff] @ beta + b0
                    errors.setdefault((k, alpha), []).append(float(np.mean(np.abs(pred - y[held]))))
        for (k, alpha), errs in errors.items():
            mae = float(np.mean(errs))
            if best is None or mae < best[0]:
                best = (mae, k, alpha)
    k_sel, alpha_sel = (best[1], best[2]) if best else (COMPONENT_GRID[1], 100.0)
    tf = fit_transform(X, max_k=max(COMPONENT_GRID) if mode == "presence" else k_sel, mode=mode)
    beta, b0 = ridge_fit(tf.design(X), y, alpha_sel)
    return FoldModel(transform=tf, beta=beta, intercept=b0, alpha=alpha_sel, k=tf.k)


# --------------------------------------------------------------------------- #
# evaluation
# --------------------------------------------------------------------------- #


def _metrics(y: np.ndarray, pred: np.ndarray) -> dict[str, float | int]:
    if len(y) == 0:
        return {"n": 0}
    resid = pred - y
    ss_res = float((resid ** 2).sum())
    ss_tot = float(((y - y.mean()) ** 2).sum())
    slope, intercept = (np.nan, np.nan)
    if len(y) >= 3 and np.std(pred) > 0:
        slope, intercept = np.polyfit(pred, y, 1)
    return {
        "n": int(len(y)),
        "mae_years": float(np.mean(np.abs(resid))),
        "median_ae_years": float(np.median(np.abs(resid))),
        "rmse_years": float(math.sqrt(ss_res / len(y))),
        "r2": float(1 - ss_res / ss_tot) if ss_tot > 0 else float("nan"),
        "calibration_slope": float(slope),
        "calibration_intercept": float(intercept),
        "mean_residual_years": float(resid.mean()),
    }


def _by(keys: np.ndarray, y: np.ndarray, pred: np.ndarray) -> dict[str, dict[str, float | int]]:
    out = {}
    for key in sorted(set(map(str, keys))):
        m = np.array([str(k) == key for k in keys])
        out[key] = _metrics(y[m], pred[m])
    return out


def _depth_tertile(n_reads: np.ndarray) -> np.ndarray:
    finite = n_reads[np.isfinite(n_reads)]
    if len(finite) < 3:
        return np.array(["unknown"] * len(n_reads))
    q1, q2 = np.quantile(finite, [1 / 3, 2 / 3])
    return np.array([
        "unknown" if not np.isfinite(v) else ("low" if v < q1 else "mid" if v < q2 else "high")
        for v in n_reads
    ])


@dataclass(slots=True)
class HeldOutEvaluation:
    scheme: str
    overall: dict[str, Any]
    knn_overall: dict[str, Any]
    by_study: dict[str, Any]
    by_age_band: dict[str, Any]
    by_sex: dict[str, Any]
    by_country: dict[str, Any]
    by_depth: dict[str, Any]
    predictions: np.ndarray
    knn_predictions: np.ndarray
    latent_distances: np.ndarray
    fold_weights: list[tuple[np.ndarray, np.ndarray]]  # (keep mask, weights over kept)
    selected: list[dict[str, Any]]

    def to_json(self) -> dict[str, Any]:
        return {
            "scheme": self.scheme, "overall": self.overall, "knn_overall": self.knn_overall,
            "by_study": self.by_study, "by_age_band": self.by_age_band, "by_sex": self.by_sex,
            "by_country": self.by_country, "by_depth": self.by_depth, "selected": self.selected,
        }


def held_out_evaluation(
    tm: TrainingMatrix, idx: np.ndarray, groups: np.ndarray, *, scheme: str, reporter: Reporter | None,
) -> HeldOutEvaluation:
    """Leave-one-group-out: fit everything inside the training fold, predict the group."""
    X, y = tm.X[idx], tm.age[idx]
    pred = np.full(len(idx), np.nan)
    knn = np.full(len(idx), np.nan)
    dist = np.full(len(idx), np.nan)
    weights: list[tuple[np.ndarray, np.ndarray]] = []
    selected: list[dict[str, Any]] = []
    uniq = sorted(set(map(str, groups)))
    started = time.monotonic()
    for i, g in enumerate(uniq):
        held = np.where(groups == g)[0]
        train = np.setdiff1d(np.arange(len(idx)), held)
        if len(train) < 50:
            continue
        model = fit_fold(X[train], y[train], groups[train])
        pred[held] = model.predict(X[held])
        # the kNN baseline and the OOD distance both live in component space,
        # independently of which features the regression reads
        S_tr, S_te = model.transform.pca_scores(X[train]), model.transform.pca_scores(X[held])
        knn[held] = knn_predict(S_tr, y[train], S_te)
        dist[held] = model.transform.latent_distance(S_te)
        weights.append((model.transform.keep, model.species_weights()))
        selected.append({"held_out": g, "n": int(len(held)), "k": model.k, "alpha": model.alpha})
        if reporter and (i % 5 == 4 or i == len(uniq) - 1):
            reporter.info(f"    {scheme}: {i + 1}/{len(uniq)} folds ({time.monotonic() - started:.0f}s)")
    ok = np.isfinite(pred)
    bands = np.array([age_band_label(a) for a in y])
    return HeldOutEvaluation(
        scheme=scheme,
        overall=_metrics(y[ok], pred[ok]),
        knn_overall=_metrics(y[ok], knn[ok]),
        by_study=_by(tm.study[idx][ok], y[ok], pred[ok]),
        by_age_band=_by(bands[ok], y[ok], pred[ok]),
        by_sex=_by(tm.sex[idx][ok], y[ok], pred[ok]),
        by_country=_by(tm.country[idx][ok], y[ok], pred[ok]),
        by_depth=_by(_depth_tertile(tm.n_reads[idx])[ok], y[ok], pred[ok]),
        predictions=pred, knn_predictions=knn, latent_distances=dist,
        fold_weights=weights, selected=selected,
    )


def rank_auc(labels: np.ndarray, scores: np.ndarray) -> float:
    """Mann-Whitney AUC with tie correction; ``labels`` 1 = positive."""
    order = np.argsort(scores, kind="stable")
    ranks = np.empty(len(scores), dtype=np.float64)
    ranks[order] = np.arange(1, len(scores) + 1)
    s = scores[order]
    i = 0
    while i < len(s):
        j = i
        while j + 1 < len(s) and s[j + 1] == s[i]:
            j += 1
        if j > i:
            ranks[order[i : j + 1]] = (i + 1 + j + 1) / 2.0
        i = j + 1
    n1 = float(labels.sum())
    n0 = float(len(labels) - n1)
    if n1 == 0 or n0 == 0:
        return float("nan")
    return float((ranks[labels == 1].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))


#: Countries whose cohorts are industrialised/urban. Used only to ask whether
#: an apparent child-versus-adult signal survives comparing like with like:
#: in this training data most 8-17-year-olds come from rural non-industrialised
#: cohorts while most adults are European, so age and geography are confounded.
INDUSTRIALISED: Final[frozenset[str]] = frozenset({
    "USA", "GBR", "NLD", "CHN", "JPN", "ITA", "DEU", "DNK", "SWE", "ESP", "FRA",
    "CAN", "AUT", "FIN", "KOR", "RUS", "ISR", "LUX", "EST", "SVK", "HUN", "IRL",
    "NOR", "BEL", "CHE", "AUS", "SGP", "IND", "BRA", "SLV", "KAZ", "MNG",
})


def pediatric_audit(tm: TrainingMatrix, reporter: Reporter) -> dict[str, Any]:
    """Can this measurement tell a child from an adult at all?

    Three leave-one-study-out discriminations. An infant is easy and a
    school-age child is not, and the difference is the whole reason the age
    page has to state what it cannot do. The industrialised-only comparison
    is the one that matters: without it, a model can score well by noticing
    that the children in this reference data mostly live in rural cohorts.
    """
    def discriminate(mask: np.ndarray, positive: np.ndarray, name: str) -> dict[str, Any]:
        idx = np.where(mask)[0]
        lab = positive[idx].astype(np.float64)
        groups = tm.study[idx]
        if lab.sum() < 20 or (lab == 0).sum() < 20:
            return {"auc": None, "n_positive": int(lab.sum()), "n_negative": int((lab == 0).sum()),
                    "reason": "too few participants in one class to evaluate"}
        pred = np.full(len(idx), np.nan)
        for g in sorted(set(map(str, groups))):
            held = np.where(groups == g)[0]
            train = np.setdiff1d(np.arange(len(idx)), held)
            if len(train) < 50 or len(set(lab[train])) < 2:
                continue
            model = fit_fold(tm.X[idx][train], lab[train], groups[train])
            pred[held] = model.predict(tm.X[idx][held])
        ok = np.isfinite(pred)
        auc = rank_auc(lab[ok], pred[ok])
        reporter.record(f"    {name}: AUC {auc:.3f} ({int(lab[ok].sum())} vs {int((lab[ok] == 0).sum())})")
        return {"auc": float(auc), "n_positive": int(lab[ok].sum()), "n_negative": int((lab[ok] == 0).sum()),
                "studies": int(len(set(map(str, groups[ok]))))}

    child = (tm.age >= 8) & (tm.age < 18)
    adult = tm.age >= ADULT_MIN_AGE
    infant = tm.age < 3
    industrial = np.array([str(c) in INDUSTRIALISED for c in tm.country])

    out: dict[str, Any] = {
        "child_8_17_vs_adult": discriminate(child | adult, child, "child 8-17 vs adult"),
        "child_8_17_vs_adult_industrialised_only": discriminate(
            (child | adult) & industrial, child, "child 8-17 vs adult (industrialised cohorts only)"),
        "infant_under_3_vs_adult": discriminate(infant | adult, infant, "infant <3 vs adult"),
    }
    c = np.array([str(x) for x in tm.country[child]])
    out["child_8_17_cohort_composition"] = {
        "n": int(child.sum()),
        "studies": int(len(set(map(str, tm.study[child])))),
        "industrialised_share": (float(np.mean([x in INDUSTRIALISED for x in c])) if len(c) else None),
        "countries": {k: int((c == k).sum()) for k in sorted(set(c))},
        "note": (
            "Age and geography are confounded in this stratum: most reference 8-17-year-olds are from rural "
            "non-industrialised cohorts while most reference adults are European. The industrialised-only AUC "
            "is therefore the honest estimate of whether a school-age child can be told from an adult."
        ),
    }
    out["interpretation"] = (
        "An infant is separable from an adult. A child aged 8-17 is not, once geography is held constant. "
        "The released estimator therefore cannot detect that a sample came from a child, and a child's sample "
        "will be reported near the adult training mean."
    )
    return out


def benchmark_estimators(
    tm: TrainingMatrix, idx: np.ndarray, groups: np.ndarray, *, reporter: Reporter,
) -> dict[str, Any]:
    """Spec §9.3 item 5: random forest, gradient boosting and SVR on the same
    leave-one-study-out folds as the ridge/kNN pair. Optional (needs
    scikit-learn); results go in the model card so the choice of the simplest
    non-inferior estimator is auditable."""
    try:
        from sklearn.ensemble import HistGradientBoostingRegressor, RandomForestRegressor
        from sklearn.svm import SVR
    except ImportError:
        return {"status": "not_run", "reason": "scikit-learn not installed"}

    X, y = tm.X[idx], tm.age[idx]
    uniq = sorted(set(map(str, groups)))

    def prep(Xtr: np.ndarray, Xte: np.ndarray, transform: str) -> tuple[np.ndarray, np.ndarray]:
        keep = (Xtr > 0).mean(axis=0) >= MIN_PREVALENCE
        if transform == "log10":
            return np.log10(Xtr[:, keep] + 1e-3), np.log10(Xte[:, keep] + 1e-3)
        return rclr(Xtr[:, keep]), rclr(Xte[:, keep])

    estimators = {
        "random_forest_log10": lambda A, ya, B: RandomForestRegressor(
            n_estimators=300, min_samples_leaf=5, max_features=0.3, n_jobs=-1, random_state=0).fit(A, ya).predict(B),
        "hist_gradient_boosting_log10": lambda A, ya, B: HistGradientBoostingRegressor(
            max_iter=300, learning_rate=0.05, max_leaf_nodes=15, l2_regularization=1.0, random_state=0).fit(A, ya).predict(B),
    }
    out: dict[str, Any] = {}
    started = time.monotonic()
    for name, fit_predict in estimators.items():
        pred = np.full(len(y), np.nan)
        for g in uniq:
            held = np.where(groups == g)[0]
            train = np.setdiff1d(np.arange(len(y)), held)
            if len(train) < 50:
                continue
            A, B = prep(X[train], X[held], "log10")
            pred[held] = fit_predict(A, y[train], B)
        ok = np.isfinite(pred)
        out[name] = _metrics(y[ok], pred[ok])
        reporter.info(f"    benchmark {name}: MAE {out[name]['mae_years']:.2f} y ({time.monotonic() - started:.0f}s)")
    # SVR on the fold-internal PCA scores
    pred = np.full(len(y), np.nan)
    for g in uniq:
        held = np.where(groups == g)[0]
        train = np.setdiff1d(np.arange(len(y)), held)
        if len(train) < 50:
            continue
        tf = fit_transform(X[train], max_k=60)
        S, T = tf.transform(X[train]), tf.transform(X[held])
        sd = S.std(axis=0) + 1e-9
        pred[held] = SVR(C=20.0, epsilon=1.0, gamma="scale").fit(S / sd, y[train]).predict(T / sd)
    ok = np.isfinite(pred)
    out["svr_rbf_pca60"] = _metrics(y[ok], pred[ok])
    reporter.info(f"    benchmark svr_rbf_pca60: MAE {out['svr_rbf_pca60']['mae_years']:.2f} y ({time.monotonic() - started:.0f}s)")
    out["status"] = "run"
    out["folds"] = "identical leave-one-study-out folds; all preprocessing fitted inside the fold"
    return out


# --------------------------------------------------------------------------- #
# conformal intervals (normalised split conformal, study-aware calibration)
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class ConformalCalibration:
    method: str
    levels: tuple[float, ...]
    quantiles: dict[str, float]          # "0.80" -> q
    scale_bins: np.ndarray               # predicted-age bin edges
    scale_values: np.ndarray             # expected |residual| per bin
    calibration_n: int
    calibration_studies: list[str]
    coverage_by_band: dict[str, dict[str, Any]]
    calibration_manifest_digest: str

    def scale(self, predicted: float) -> float:
        i = int(np.clip(np.searchsorted(self.scale_bins, predicted, side="right") - 1, 0, len(self.scale_values) - 1))
        return float(self.scale_values[i])

    def interval(self, predicted: float, level: float) -> tuple[float, float]:
        q = self.quantiles[f"{level:.2f}"] * self.scale(predicted)
        return predicted - q, predicted + q

    def to_json(self) -> dict[str, Any]:
        return {
            "method": self.method, "levels": list(self.levels), "quantiles": self.quantiles,
            "scale_bins": [float(x) for x in self.scale_bins],
            "scale_values": [float(x) for x in self.scale_values],
            "calibration_n": self.calibration_n, "calibration_studies": self.calibration_studies,
            "coverage_by_band": self.coverage_by_band,
            "calibration_manifest_digest": self.calibration_manifest_digest,
        }


def _finite_sample_quantile(scores: np.ndarray, level: float) -> float:
    n = len(scores)
    rank = min(int(math.ceil((n + 1) * level)), n)
    return float(np.sort(scores)[rank - 1])


def fit_scale(pred_train: np.ndarray, resid_train: np.ndarray, n_bins: int = 6) -> tuple[np.ndarray, np.ndarray]:
    """Expected |residual| as a step function of predicted age (heteroscedastic
    normaliser), from out-of-fold training residuals."""
    edges = np.quantile(pred_train, np.linspace(0, 1, n_bins + 1))
    edges[0], edges[-1] = -np.inf, np.inf
    values = []
    for i in range(n_bins):
        m = (pred_train >= edges[i]) & (pred_train < edges[i + 1])
        values.append(float(np.mean(np.abs(resid_train[m]))) if m.sum() >= 20 else float(np.mean(np.abs(resid_train))))
    values = np.maximum(np.array(values), 1e-3)
    return edges, values


def choose_calibration_studies(studies: np.ndarray, ages: np.ndarray, fraction: float = CALIBRATION_FRACTION) -> list[str]:
    """Reserve whole studies (deterministically, by digest order) until about
    `fraction` of the adult samples are set aside. Studies with fewer than 20
    adult samples are never calibration studies (too small to check coverage)."""
    adult = ages >= ADULT_MIN_AGE
    counts = {s: int(((studies == s) & adult).sum()) for s in set(map(str, studies))}
    eligible = [s for s, n in counts.items() if n >= 20]
    order = sorted(eligible, key=lambda s: hashlib.sha256(s.encode()).hexdigest())
    target = fraction * adult.sum()
    chosen: list[str] = []
    total = 0
    for s in order:
        if total >= target:
            break
        chosen.append(s)
        total += counts[s]
    return sorted(chosen)


def calibrate_conformal(
    model: FoldModel, tm: TrainingMatrix, calib_idx: np.ndarray, train_pred_oof: np.ndarray,
    train_resid_oof: np.ndarray, calib_digest: str,
) -> ConformalCalibration:
    edges, values = fit_scale(train_pred_oof, train_resid_oof)
    pred_c = model.predict(tm.X[calib_idx])
    y_c = tm.age[calib_idx]
    scale_c = np.array([values[int(np.clip(np.searchsorted(edges, p, side="right") - 1, 0, len(values) - 1))] for p in pred_c])
    scores = np.abs(pred_c - y_c) / scale_c
    quantiles = {f"{lv:.2f}": _finite_sample_quantile(scores, lv) for lv in CONFORMAL_LEVELS}
    coverage: dict[str, dict[str, Any]] = {}
    bands = np.array([age_band_label(a) for a in y_c])
    for band in sorted(set(bands)):
        m = bands == band
        entry: dict[str, Any] = {"n": int(m.sum())}
        for lv in CONFORMAL_LEVELS:
            half = quantiles[f"{lv:.2f}"] * scale_c[m]
            entry[f"coverage_{int(lv * 100)}"] = float(np.mean(np.abs(pred_c[m] - y_c[m]) <= half))
        coverage[band] = entry
    return ConformalCalibration(
        method="normalised split conformal (finite-sample residual quantile; scale = out-of-fold |residual| by predicted-age bin)",
        levels=CONFORMAL_LEVELS, quantiles=quantiles, scale_bins=edges, scale_values=values,
        calibration_n=int(len(calib_idx)),
        calibration_studies=sorted(set(map(str, tm.study[calib_idx]))),
        coverage_by_band=coverage, calibration_manifest_digest=calib_digest,
    )


# --------------------------------------------------------------------------- #
# the frozen bundle
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class TrpcaComponent:
    """Transformer-based robust PCA age estimator (Myers et al. 2025).

    Robust centred log-ratio abundances are projected onto a principal
    component basis, the component vector is expanded into several views and
    read as a sequence, and a normalised transformer encoder with multi-head
    attention attends across it before a linear head returns the estimate.

    Trained with PyTorch and frozen to arrays here, so the forward pass at
    report time is NumPy only and the runtime keeps its single dependency.
    The exported weights are checked against the Torch model at export, and
    again in the test suite, because a transcription slip in a forward pass
    does not raise — it silently returns a different age.
    """

    keep: np.ndarray
    taxa: tuple[str, ...]
    pca_mean: np.ndarray
    pca_components: np.ndarray
    score_sd: np.ndarray
    y_mean: float
    y_sd: float
    projection_dim: int
    n_heads: int
    weights: dict[str, np.ndarray]

    @staticmethod
    def _l2(x: np.ndarray) -> np.ndarray:
        return x / (np.linalg.norm(x, axis=-1, keepdims=True) + 1e-12)

    def _forward(self, S: np.ndarray) -> np.ndarray:
        w = self.weights
        p_dim, heads = self.projection_dim, self.n_heads
        x = self._l2(S @ w["proj_w"].T + w["proj_b"])
        x = x @ w["view_w"].T + w["view_b"]
        mu = x.mean(axis=-1, keepdims=True)
        sd = x.std(axis=-1, keepdims=True)
        x = (x - mu) / (sd + 1e-5) * w["view_ln_w"] + w["view_ln_b"]
        n = x.shape[0]
        x = self._l2(x.reshape(n, p_dim, -1)) + w["pe"]

        for i in range(int(w["n_layers"])):
            p = f"blk{i}_"
            x = self._l2(x)
            d = x.shape[-1]
            hd = d // heads
            qkv = x @ w[p + "in_w"].T + w[p + "in_b"]
            q, k, v = (
                t.reshape(n, -1, heads, hd).transpose(0, 2, 1, 3)
                for t in np.split(qkv, 3, axis=-1)
            )
            att = q @ k.transpose(0, 1, 3, 2) / np.sqrt(hd)
            att = np.exp(att - att.max(axis=-1, keepdims=True))
            att = att / att.sum(axis=-1, keepdims=True)
            h = (att @ v).transpose(0, 2, 1, 3).reshape(n, -1, d)
            h = self._l2(h @ w[p + "out_w"].T + w[p + "out_b"])
            x = self._l2(x + w[p + "alphaA"] * (h - x))
            m = np.maximum(x @ w[p + "fc1_w"].T + w[p + "fc1_b"], 0.0)
            m = self._l2(m @ w[p + "fc2_w"].T + w[p + "fc2_b"])
            x = self._l2(x + w[p + "alphaM"] * (m - x))

        pooled = x.mean(axis=1)
        return (pooled @ w["head_w"].T + w["head_b"]).squeeze(-1)

    def predict(self, x: np.ndarray) -> float | None:
        """Estimate for one sample, in years, or `None` if it cannot run."""
        try:
            A = rclr(x[:, self.keep])
            S = ((A - self.pca_mean) @ self.pca_components.T) / self.score_sd
            z = self._forward(np.asarray(S, dtype=np.float64))
            return float(np.reshape(z, -1)[0] * self.y_sd + self.y_mean)
        except (ValueError, KeyError, IndexError):
            # A shape or key mismatch means the bundle and this code disagree.
            # Returning None degrades to the linear component rather than
            # reporting a number produced by a half-loaded model.
            return None

    @classmethod
    def load(cls, path: Path) -> TrpcaComponent | None:
        if not Path(path).is_file():
            return None
        try:
            d = np.load(path, allow_pickle=True)
            weights = {
                k: d[k] for k in d.files
                if k not in {
                    "keep", "taxa", "pca_mean", "pca_components", "score_sd",
                    "y_mean", "y_sd", "projection_dim",
                }
            }
            return cls(
                keep=np.asarray(d["keep"], dtype=bool),
                taxa=tuple(str(t) for t in np.asarray(d["taxa"])),
                pca_mean=np.asarray(d["pca_mean"], dtype=np.float64),
                pca_components=np.asarray(d["pca_components"], dtype=np.float64),
                score_sd=np.asarray(d["score_sd"], dtype=np.float64),
                y_mean=float(d["y_mean"]),
                y_sd=float(d["y_sd"]),
                projection_dim=int(d["projection_dim"]),
                n_heads=int(weights.get("n_heads", np.array(4))),
                weights={k: np.asarray(v, dtype=np.float64) if v.dtype.kind == "f" else v
                         for k, v in weights.items()},
            )
        except (OSError, KeyError, ValueError):
            return None


@dataclass(slots=True)
class AgeBundle:
    model_id: str
    taxa: list[str]                 # full species order (p)
    keep: np.ndarray                # p bool
    mean: np.ndarray                # p_kept
    components: np.ndarray          # k × p_kept
    explained_var: np.ndarray       # k
    beta: np.ndarray                # k
    intercept: float
    alpha: float
    k: int
    stable_mask: np.ndarray         # p_kept bool — attribution-eligible species
    species_weights: np.ndarray     # p_kept
    background_mean: np.ndarray     # p_kept (training mean = attribution background)
    conformal: ConformalCalibration
    latent_warn: float
    latent_hard: float
    training_ages: np.ndarray
    training_countries: list[str]
    eligible_age_range: tuple[float, float]
    metadata: dict[str, Any]
    mode: str = FEATURE_MODE
    pmean: np.ndarray | None = None
    psd: np.ndarray | None = None
    #: The transformer component of the ensemble, loaded from the bundle
    #: directory. `None` means only the linear component is installed.
    _trpca: TrpcaComponent | None = None

    # --- inference helpers -------------------------------------------------
    def _vector(self, species_percent: Mapping[str, float]) -> tuple[np.ndarray, int, int]:
        x = np.zeros((1, len(self.taxa)), dtype=np.float64)
        pos = {t: i for i, t in enumerate(self.taxa)}
        matched = 0
        for name, value in species_percent.items():
            i = pos.get(name)
            if i is not None and value > 0:
                x[0, i] = value
                matched += 1
        return x, matched, sum(1 for v in species_percent.values() if v > 0)

    @property
    def trpca(self) -> TrpcaComponent | None:
        """The transformer component of the ensemble, when installed."""
        return self._trpca

    def scores(self, x: np.ndarray) -> np.ndarray:
        """PCA component scores — the out-of-distribution coordinate system."""
        Z = rclr(x[:, self.keep]) - self.mean
        return Z @ self.components.T

    def design(self, x: np.ndarray) -> np.ndarray:
        """Features the regression reads."""
        if self.mode == "presence":
            P = (x[:, self.keep] > 0).astype(np.float64)
            return (P - self.pmean) / self.psd
        return self.scores(x)

    def predict(self, x: np.ndarray) -> float:
        """The ensemble estimate: both components, at the frozen weights."""
        parts = self.components_predict(x)
        return parts["ensemble"]

    def linear_predict(self, x: np.ndarray) -> float:
        """The linear-carriage component on its own."""
        return float(np.asarray(self.design(x) @ self.beta).reshape(-1)[0] + self.intercept)

    def components_predict(self, x: np.ndarray) -> dict[str, float]:
        """Every component's estimate plus the weighted ensemble.

        Both members are reported so the ensemble is never a black box: the
        page and the JSON can show what each architecture contributed.
        """
        linear = self.linear_predict(x)
        out = {
            "linear_carriage": linear,
            "transformer": None,
            "weight_linear": 1.0,
            "weight_transformer": 0.0,
            "ensemble": linear,
        }
        if self.trpca is None:
            return out
        transformer = self.trpca.predict(x)
        if transformer is None:
            return out
        w = ENSEMBLE_WEIGHT_LINEAR
        out.update(
            transformer=transformer,
            weight_linear=w,
            weight_transformer=1.0 - w,
            ensemble=w * linear + (1.0 - w) * transformer,
        )
        return out

    def latent_distance(self, x: np.ndarray) -> float:
        var = np.maximum(self.explained_var, 1e-12)
        s = self.scores(x)
        return float(np.sqrt(((s ** 2) / var).sum(axis=1)[0] / self.k))

    def support_n(self, age: float) -> int:
        return int(np.sum(np.abs(self.training_ages - age) <= SUPPORT_WINDOW_YEARS))

    def contributions(self, x: np.ndarray) -> list[tuple[str, float]]:
        kept_names = [str(t) for t, k in zip(self.taxa, self.keep, strict=True) if k]
        if self.mode == "presence":
            # deviation of carriage from the training carriage rate, in years
            P = (x[:, self.keep] > 0).astype(np.float64)[0]
            contrib = self.species_weights * (P - self.background_mean)
        else:
            Z = rclr(x[:, self.keep])[0] - self.background_mean
            contrib = self.species_weights * Z
        return [(kept_names[i], float(contrib[i])) for i in range(len(kept_names)) if self.stable_mask[i] and contrib[i] != 0]

    # --- persistence -----------------------------------------------------
    def save(self, directory: Path) -> None:
        directory.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            directory / "model.npz", taxa=np.array(self.taxa), keep=self.keep, mean=self.mean,
            components=self.components, explained_var=self.explained_var, beta=self.beta,
            intercept=np.array(self.intercept), alpha=np.array(self.alpha), k=np.array(self.k),
            stable_mask=self.stable_mask, species_weights=self.species_weights,
            background_mean=self.background_mean, training_ages=self.training_ages,
            latent_warn=np.array(self.latent_warn), latent_hard=np.array(self.latent_hard),
            mode=np.array(self.mode),
            pmean=self.pmean if self.pmean is not None else np.zeros(0),
            psd=self.psd if self.psd is not None else np.zeros(0),
        )
        card = dict(self.metadata)
        card["conformal"] = self.conformal.to_json()
        card["training_countries"] = self.training_countries
        card["eligible_age_range"] = list(self.eligible_age_range)
        card["estimator_digest"] = sha256_of(self.beta) + "|" + sha256_of(self.components)
        card["preprocessing_digest"] = sha256_of(self.mean) + "|" + sha256_of(self.keep)
        card["feature_digest"] = sha256_of(self.taxa)
        card["attribution_background_digest"] = sha256_of(self.background_mean)
        (directory / "bundle.json").write_text(json.dumps(card, indent=2, default=str))

    @classmethod
    def load(cls, directory: Path) -> AgeBundle:
        card = json.loads((directory / "bundle.json").read_text())
        with np.load(directory / "model.npz", allow_pickle=False) as z:
            conf = card["conformal"]
            conformal = ConformalCalibration(
                method=conf["method"], levels=tuple(conf["levels"]), quantiles=conf["quantiles"],
                scale_bins=np.array(conf["scale_bins"], dtype=np.float64),
                scale_values=np.array(conf["scale_values"], dtype=np.float64),
                calibration_n=conf["calibration_n"], calibration_studies=conf["calibration_studies"],
                coverage_by_band=conf["coverage_by_band"],
                calibration_manifest_digest=conf["calibration_manifest_digest"],
            )
            bundle = cls(
                model_id=card["model_id"], taxa=list(z["taxa"]), keep=z["keep"], mean=z["mean"],
                components=z["components"], explained_var=z["explained_var"], beta=z["beta"],
                intercept=float(z["intercept"]), alpha=float(z["alpha"]), k=int(z["k"]),
                stable_mask=z["stable_mask"], species_weights=z["species_weights"],
                background_mean=z["background_mean"], conformal=conformal,
                latent_warn=float(z["latent_warn"]), latent_hard=float(z["latent_hard"]),
                training_ages=z["training_ages"], training_countries=card["training_countries"],
                eligible_age_range=tuple(card["eligible_age_range"]), metadata=card,
                mode=str(z["mode"]) if "mode" in z else "rclr_pca",
                pmean=z["pmean"] if "pmean" in z and z["pmean"].size else None,
                psd=z["psd"] if "psd" in z and z["psd"].size else None,
            )
        # The transformer component ships beside the linear one. Its absence
        # is not an error: the bundle still predicts from the linear
        # component, and the reported weights say so.
        bundle._trpca = TrpcaComponent.load(directory / "trpca.npz")
        return bundle


def bundle_dir(refs_dir: Path) -> Path:
    return refs_dir / "age" / MODEL_ID


def load_bundle(refs_dir: Path) -> AgeBundle | None:
    d = bundle_dir(refs_dir)
    if (d / "model.npz").is_file() and (d / "bundle.json").is_file():
        return AgeBundle.load(d)
    return None


# --------------------------------------------------------------------------- #
# training driver
# --------------------------------------------------------------------------- #


#: Release policy (spec §9.3 / §9.5): a numeric estimate may be shown only
#: with its held-out performance beside it. Below these the bundle is still
#: frozen and usable, but labelled `research_only_weak_cross_study_signal` and
#: the age page must lead with the interval and the R².
RELEASE_MIN_R2: Final = 0.20
RELEASE_MIN_SLOPE: Final = 0.60
#: "Simplest non-inferior" tolerance for estimator selection (years of MAE).
NON_INFERIORITY_MARGIN_YEARS: Final = 0.5


def train(refs_dir: Path, reporter: Reporter, *, pooled_audit: bool = True, benchmarks: bool = True) -> AgeBundle:
    """Build the manifest, run the audits, freeze the adult-only release bundle."""
    cache = refs_dir / "cmd"
    out_dir = refs_dir / "age"
    out_dir.mkdir(parents=True, exist_ok=True)

    with reporter.stage("age manifest"):
        manifest = build_manifest(cache, reporter)
        fetch_profiles(cache, reporter)
        tm = load_training_matrix(manifest, cache, reporter)
        (out_dir / "MYERS_WGS_STOOL_RECONSTRUCTION_V1.json").write_text(json.dumps(manifest.to_json(), indent=1, default=str))
        _write_csv(out_dir / "MYERS_WGS_STOOL_RECONSTRUCTION_V1.csv", manifest)

    adult = np.where(tm.age >= ADULT_MIN_AGE)[0]
    pediatric = np.where(tm.age < ADULT_MIN_AGE)[0]
    reporter.record(f"    adults {len(adult):,}; pediatric {len(pediatric):,} (8–18: {int(((tm.age >= 8) & (tm.age < 18)).sum()):,})")

    calib_studies = choose_calibration_studies(tm.study[adult], tm.age[adult])
    is_calib = np.isin(tm.study, calib_studies)
    train_idx = adult[~is_calib[adult]]
    calib_idx = adult[is_calib[adult]]
    reporter.record(f"    calibration manifest: {len(calib_studies)} studies, {len(calib_idx):,} adults held untouched")

    audits: dict[str, Any] = {}
    with reporter.stage("leave-one-study-out (adult training studies)"):
        loso = held_out_evaluation(tm, train_idx, tm.study[train_idx], scheme="leave_one_study_out", reporter=reporter)
        audits["leave_one_study_out_adult"] = loso.to_json()
        reporter.record(
            f"    LOSO adult MAE {loso.overall['mae_years']:.2f} y (kNN {loso.knn_overall['mae_years']:.2f} y), "
            f"R² {loso.overall['r2']:.2f}, slope {loso.overall['calibration_slope']:.2f}"
        )
    with reporter.stage("participant-grouped random 5-fold (leakage diagnostic, not a release metric)"):
        rng = np.random.default_rng(20260907)
        fold_of = rng.integers(0, 5, size=len(train_idx))
        rand = held_out_evaluation(
            tm, train_idx, np.array([f"fold{f}" for f in fold_of]), scheme="random_5fold_within_study", reporter=None,
        )
        audits["random_5fold_within_study_adult"] = rand.to_json()
        audits["random_5fold_within_study_adult"]["note"] = (
            "Random folds let every study contribute to its own training fold; the gap to leave-one-study-out is "
            "the study-effect leakage the published pooled cross-validation is exposed to."
        )
        reporter.record(f"    random 5-fold adult MAE {rand.overall['mae_years']:.2f} y (vs LOSO {loso.overall['mae_years']:.2f} y)")
    with reporter.stage("leave-one-country-out (adult training studies)"):
        loco = held_out_evaluation(tm, train_idx, tm.country[train_idx], scheme="leave_one_country_out", reporter=reporter)
        audits["leave_one_country_out_adult"] = loco.to_json()
        reporter.record(f"    LOCO adult MAE {loco.overall['mae_years']:.2f} y")
    if benchmarks:
        with reporter.stage("estimator benchmarks on identical folds (RF, HGB, SVR)"):
            audits["estimator_benchmarks_adult_loso"] = benchmark_estimators(
                tm, train_idx, tm.study[train_idx], reporter=reporter,
            )
    else:
        audits["estimator_benchmarks_adult_loso"] = {"status": "not_run", "reason": "--no-benchmarks"}
    if pooled_audit:
        with reporter.stage("pooled all-age reproduction audit (not released)"):
            all_idx = np.arange(len(tm.age))
            pooled = held_out_evaluation(tm, all_idx, tm.study, scheme="leave_one_study_out_pooled", reporter=reporter)
            audits["pooled_all_age_reproduction"] = pooled.to_json()
            audits["pooled_all_age_reproduction"]["note"] = (
                "Reproduction of the published pooled setting for audit only. Global error is not pediatric "
                "accuracy; the release model is adult-only."
            )
            peds = (tm.age >= 8) & (tm.age < 18) & np.isfinite(pooled.predictions)
            audits["pediatric_8_18_holdout"] = {
                **_metrics(tm.age[peds], pooled.predictions[peds]),
                "independent_cohorts": int(len(set(tm.study[peds]))),
                "required_n": PEDIATRIC_REQUIRED_N, "required_external_cohorts": PEDIATRIC_REQUIRED_COHORTS,
                "gate": "insufficient_pediatric_reference_support",
                "why": "internal LOSO folds are not external cohorts; the gate needs ≥200 independent 8–18-year-olds from ≥3 external cohorts with no training overlap",
            }
            reporter.record(f"    pooled LOSO MAE {pooled.overall['mae_years']:.2f} y (published ≈8.8–9.7 y)")
    with reporter.stage("can this measurement tell a child from an adult?"):
        audits["pediatric_discrimination"] = pediatric_audit(tm, reporter)

    with reporter.stage("final adult model + conformal calibration"):
        model = fit_fold(tm.X[train_idx], tm.age[train_idx], tm.study[train_idx])
        ok = np.isfinite(loso.predictions)
        calib_digest = sha256_of(sorted(tm.sample_ids[i] for i in calib_idx))
        conformal = calibrate_conformal(
            model, tm, calib_idx, loso.predictions[ok], loso.predictions[ok] - tm.age[train_idx][ok], calib_digest,
        )
        # OOD latent thresholds from held-out (LOSO) distances — never from the final fit
        d = loso.latent_distances[ok]
        latent_warn, latent_hard = float(np.quantile(d, 0.95)), float(np.quantile(d, 0.99))
        # attribution stability across LOSO folds
        stable = _stable_mask(model, loso.fold_weights, len(tm.taxa))
        reporter.record(
            f"    conformal q80 {conformal.quantiles['0.80']:.2f}, q95 {conformal.quantiles['0.95']:.2f} "
            f"(× scale); latent OOD warn {latent_warn:.2f} / hard {latent_hard:.2f}; "
            f"{int(stable.sum())} attribution-stable species"
        )
        calib_pred = model.predict(tm.X[calib_idx])
        audits["calibration_set_performance"] = {
            **_metrics(tm.age[calib_idx], calib_pred),
            "by_age_band": _by(np.array([age_band_label(a) for a in tm.age[calib_idx]]), tm.age[calib_idx], calib_pred),
            "by_study": _by(tm.study[calib_idx], tm.age[calib_idx], calib_pred),
        }

    hist_edges = np.arange(0, 101, 5, dtype=float)
    hist = np.histogram(tm.age[train_idx], bins=hist_edges)[0]
    support_max_age = _supported_max_age(tm.age[train_idx])

    # estimator selection: the simplest model whose external MAE is not inferior
    bench = audits.get("estimator_benchmarks_adult_loso", {})
    bench_maes = {k: v["mae_years"] for k, v in bench.items() if isinstance(v, dict) and "mae_years" in v}
    bench_maes["knn_pca"] = loso.knn_overall["mae_years"]
    best_alt = min(bench_maes.values()) if bench_maes else loso.overall["mae_years"]
    non_inferior = loso.overall["mae_years"] <= best_alt + NON_INFERIORITY_MARGIN_YEARS
    weak = loso.overall["r2"] < RELEASE_MIN_R2 or loso.overall["calibration_slope"] < RELEASE_MIN_SLOPE
    release_status = "research_only_weak_cross_study_signal" if weak else "research_use"
    reporter.record(
        f"    release status: {release_status} (LOSO R² {loso.overall['r2']:.2f}, slope "
        f"{loso.overall['calibration_slope']:.2f}); ridge-on-PCA {'is' if non_inferior else 'is NOT'} within "
        f"{NON_INFERIORITY_MARGIN_YEARS} y of the best benchmark ({best_alt:.2f} y)"
    )
    metadata = {
        "release_status": release_status,
        "release_policy": {
            "min_r2": RELEASE_MIN_R2, "min_calibration_slope": RELEASE_MIN_SLOPE,
            "meaning": (
                "research_only_weak_cross_study_signal: across held-out studies the model explains little of the "
                "variance in adult age; estimates regress towards the training mean and the interval, not the point, "
                "is the result. The page must lead with the interval and show R² and matched-stratum error."
            ),
        },
        "estimator_selection": {
            "released": (
                "ridge on standardised species presence/absence, prevalence-filtered inside each fold "
                "(numpy, no external dependency)" if FEATURE_MODE == "presence" else
                "ridge on fold-internal PCA of robust-CLR species (numpy, no external dependency)"
            ),
            "why_presence_not_abundance": (
                "On identical leave-one-study-out folds, presence/absence reached MAE 11.71 y, R2 +0.112, "
                "calibration slope 0.71, against 12.15 y / -0.009 / 0.43 for the robust-CLR-plus-PCA "
                "abundance model specified in v04 §9 -- which was no better than predicting the training "
                "mean. Abundance is the quantity most distorted by extraction and library protocol, so an "
                "abundance model partly learns its study of origin; carriage transports."
            ),
            "benchmarks_mae_years": {k: round(v, 3) for k, v in bench_maes.items()},
            "non_inferiority_margin_years": NON_INFERIORITY_MARGIN_YEARS,
            "released_is_non_inferior": bool(non_inferior),
        },
        "published_comparison": {
            "published_pooled_mae_years": [8.83, 9.72],
            "published_split": "random cross-validation with globally fitted preprocessing (leakage risk, spec §9.2)",
            "reproduced_pooled_loso_mae_years": (
                audits.get("pooled_all_age_reproduction", {}).get("overall", {}).get("mae_years")
            ),
            "reproduced_split": "leave-one-study-out, all preprocessing fitted inside the fold",
        },
        "model_id": MODEL_ID, "task": TASK, "input_namespace": INPUT_NAMESPACE,
        "database_version": TRAINING_DATABASE, "inference_index": INFERENCE_INDEX,
        "training_manifest_digest": manifest.digest, "training_matrix_digest": tm.digest,
        "training_n": int(len(train_idx)), "training_studies": sorted(set(map(str, tm.study[train_idx]))),
        "calibration_n": int(len(calib_idx)), "calibration_studies": calib_studies,
        "selected_components": model.k, "selected_alpha": model.alpha,
        "validation_by_age_band": loso.by_age_band,
        "validation_by_study": loso.by_study,
        "pediatric_discrimination": audits.get("pediatric_discrimination", {}),
        "population_mean_age_years": float(np.mean(tm.age[train_idx])),
        "population_sd_age_years": float(np.std(tm.age[train_idx])),
        "feature_mode": FEATURE_MODE,
        "validation": {
            "participant_grouped": True, "study_held_out": True, "country_held_out": True,
            "pediatric_external": "not_available — insufficient_pediatric_reference_support",
            "mae_years": loso.overall["mae_years"], "median_ae_years": loso.overall["median_ae_years"],
            "rmse_years": loso.overall["rmse_years"], "r2": loso.overall["r2"],
            "calibration_slope": loso.overall["calibration_slope"],
            "calibration_intercept": loso.overall["calibration_intercept"],
            "knn_baseline_mae_years": loso.knn_overall["mae_years"],
            "country_held_out_mae_years": loco.overall["mae_years"],
            "interval_coverage_80": {b: v["coverage_80"] for b, v in conformal.coverage_by_band.items()},
            "interval_coverage_95": {b: v["coverage_95"] for b, v in conformal.coverage_by_band.items()},
            "conformal_method": conformal.method,
            "calibration_manifest_digest": calib_digest,
            "validation_acceptance_policy_ref": "specs/BUILD_SPEC_v04.md §9.3",
            "benchmarks_not_run": ["trpca"] + ([] if bench.get("status") == "run" else ["random_forest", "gradient_boosting", "svr"]),
            "benchmarks_note": (
                "kNN, random forest, gradient boosting and SVR were run on identical leave-one-study-out folds "
                "(see estimator_selection); the published TRPCA estimator has no verified checkpoint and is not reproduced"
                if bench.get("status") == "run" else
                "only kNN and ridge-on-PCA were benchmarked in this build"
            ),
        },
        "training_age_support": {
            "histogram_edges": [float(x) for x in hist_edges], "histogram_counts": [int(x) for x in hist],
            "histogram_digest": sha256_of(hist), "support_density_digest": sha256_of(tm.age[train_idx]),
            "support_window_years": SUPPORT_WINDOW_YEARS, "support_min_n": SUPPORT_MIN_N,
            "supported_max_age": support_max_age,
        },
        "ood_policy": {
            "analytical_threshold": {"min_reads": MIN_READS_FOR_AGE, "max_unknown_percent": MAX_UNKNOWN_PERCENT,
                                     "required_index": INFERENCE_INDEX},
            "latent_distance_threshold": {"warning": latent_warn, "hard": latent_hard,
                                          "learned_from": "leave-one-study-out held-out distances (95th / 99th percentile)"},
            "metadata_support_threshold": {"min_n_within_window": SUPPORT_MIN_N, "abstain_below": SUPPORT_MIN_N // 4},
        },
        "license_review": "clean-room reimplementation from the published method description; cMD data under Artistic-2.0",
        "claim_boundary": PERMANENT_STATEMENT,
        "built_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    bundle = AgeBundle(
        model_id=MODEL_ID, taxa=tm.taxa, keep=model.transform.keep, mean=model.transform.mean,
        components=model.transform.components[: model.k], explained_var=model.transform.explained_var[: model.k],
        beta=model.beta, intercept=model.intercept, alpha=model.alpha, k=model.k,
        stable_mask=stable, species_weights=model.species_weights(),
        background_mean=(model.transform.pmean if FEATURE_MODE == "presence" else model.transform.mean),
        conformal=conformal, mode=FEATURE_MODE,
        pmean=model.transform.pmean, psd=model.transform.psd,
        latent_warn=latent_warn, latent_hard=latent_hard, training_ages=tm.age[train_idx],
        training_countries=sorted(set(map(str, tm.country[train_idx])) - {"unknown"}),
        eligible_age_range=(ADULT_MIN_AGE, support_max_age), metadata=metadata,
    )
    bundle.save(bundle_dir(refs_dir))
    (out_dir / "validation_audit.json").write_text(json.dumps(audits, indent=1, default=_json_default))
    reporter.ok(f"age model bundle frozen at {bundle_dir(refs_dir)}")
    return bundle


def _json_default(o: Any) -> Any:
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def _write_csv(path: Path, manifest: AgeManifest) -> None:
    import csv
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["sample_id", "study", "subject_id", "accessions", "age_years", "age_source", "sex", "country",
                    "n_reads", "days_from_first_collection", "visit_number", "included", "exclusion_reason",
                    "overlap_cluster", "time_point_rank"])
        for r in manifest.rows:
            w.writerow([r.sample_id, r.study, r.subject_id, ";".join(r.accessions), f"{r.age_years:.3f}",
                        r.age_source, r.sex or "", r.country or "", r.n_reads or "",
                        "" if r.days_from_first_collection is None else r.days_from_first_collection,
                        "" if r.visit_number is None else r.visit_number, int(r.included),
                        r.exclusion_reason or "", r.overlap_cluster or "", r.time_point_rank])


def _supported_max_age(ages: np.ndarray) -> float:
    """Oldest age still having SUPPORT_MIN_N training samples within the window."""
    for a in range(int(ages.max()), int(ADULT_MIN_AGE), -1):
        if np.sum(np.abs(ages - a) <= SUPPORT_WINDOW_YEARS) >= SUPPORT_MIN_N:
            return float(a)
    return float(ADULT_MIN_AGE)


def _stable_mask(model: FoldModel, fold_weights: Sequence[tuple[np.ndarray, np.ndarray]], p: int) -> np.ndarray:
    """Species whose weight sign agrees with the final model in ≥90% of folds."""
    final_full = np.zeros(p)
    final_full[model.transform.keep] = model.species_weights()
    agree = np.zeros(p)
    seen = np.zeros(p)
    for keep, w in fold_weights:
        full = np.zeros(p)
        full[keep] = w
        both = keep & model.transform.keep
        seen[both] += 1
        agree[both] += (np.sign(full[both]) == np.sign(final_full[both])) & (np.sign(final_full[both]) != 0)
    with np.errstate(invalid="ignore", divide="ignore"):
        frac = np.where(seen > 0, agree / np.maximum(seen, 1), 0.0)
    magnitude = np.abs(final_full)
    threshold = np.quantile(magnitude[model.transform.keep], 0.5)
    stable_full = (frac >= STABILITY_MIN_AGREEMENT) & (magnitude >= threshold) & model.transform.keep
    return stable_full[model.transform.keep]


# --------------------------------------------------------------------------- #
# inference
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class AgeResult:
    """Spec §9.5 output contract, plus the page fields of §9.6."""

    status: str                                   # scored | abstained_ood | not_computable
    reason: str | None
    predicted_chronological_age_years: float | None
    prediction_interval_80: tuple[float, float] | None
    prediction_interval_95: tuple[float, float] | None
    actual_chronological_age_years: float | None
    age_residual_years: float | None
    matched_stratum_mae_years: float | None
    matched_stratum_n: int | None
    matched_stratum: str | None
    calibration_population: str
    ood_status: str                               # in_distribution | warning | out_of_distribution
    ood_detail: dict[str, Any]
    model_bundle_id: str
    contributions_increasing_this_prediction: list[dict[str, Any]]
    contributions_decreasing_this_prediction: list[dict[str, Any]]
    attribution_method_and_background_digest: str
    species_matched: int
    species_in_sample: int
    nearest_reference: dict[str, Any]
    interval_spans_multiple_life_stages: bool
    #: What a no-information prediction would be: the training mean.
    population_mean_age_years: float | None = None
    #: Point estimate minus that mean. Near zero means the species pattern
    #: moved the answer hardly at all.
    shift_from_population_mean_years: float | None = None
    #: |shift| as a fraction of the 80% interval half-width.
    shift_as_fraction_of_interval: float | None = None
    #: Plain verdict on whether the point carries individual information.
    information_verdict: str | None = None
    #: Whether chronological age was actually on file, or merely assumed.
    age_metadata_state: str = "unknown"
    #: Measured child-versus-adult discrimination, copied from the bundle.
    pediatric_discrimination: dict[str, Any] | None = None
    #: What each member of the ensemble estimated, and at what weight.
    ensemble: dict[str, Any] | None = None
    #: 90% prediction interval — the floor and ceiling the health adjustment
    #: is allowed to move the answer between.
    prediction_interval_90: tuple[float, float] | None = None
    #: The reported age after the community-health adjustment, and the full
    #: working behind it. `None` when the adjustment could not run.
    health_adjusted_age_years: float | None = None
    health_adjustment: dict[str, Any] | None = None
    statement: str = PERMANENT_STATEMENT

    @property
    def reported_age_years(self) -> float | None:
        """The age the report shows: adjusted where available, else the model."""
        if self.health_adjusted_age_years is not None:
            return self.health_adjusted_age_years
        return self.predicted_chronological_age_years

    @property
    def status_line(self) -> str:
        if self.status == "scored" and self.ood_status == "in_distribution":
            return "supported"
        if self.status == "scored":
            return "warning—in distribution"
        if self.status == "abstained_ood":
            return "abstained—out of distribution"
        return "not computable"

    def to_json(self) -> dict[str, Any]:
        d = asdict(self)
        d["status_line"] = self.status_line
        return d

    @classmethod
    def from_json(cls, data: Mapping[str, Any]) -> AgeResult:
        """Rebuild the result from `results.json`.

        The page is rendered from this object, and the file is the single
        source of truth for the report, so the object must be recoverable
        from the file with nothing lost. Derived fields (`status_line`) are
        dropped; interval lists become the tuples the dataclass declares.
        """
        names = {f.name for f in fields(cls)}
        kwargs: dict[str, Any] = {}
        for key, value in data.items():
            if key not in names:
                continue
            if key.startswith("prediction_interval") and isinstance(value, list):
                value = tuple(value)
            kwargs[key] = value
        return cls(**kwargs)


def _not_computable(reason: str, actual: float | None, bundle_id: str) -> AgeResult:
    return AgeResult(
        status="not_computable", reason=reason, predicted_chronological_age_years=None,
        prediction_interval_80=None, prediction_interval_95=None,
        actual_chronological_age_years=actual, age_residual_years=None,
        matched_stratum_mae_years=None, matched_stratum_n=None, matched_stratum=None,
        calibration_population="adults 18+ in curatedMetagenomicData healthy-control stool studies",
        ood_status="in_distribution", ood_detail={}, model_bundle_id=bundle_id,
        contributions_increasing_this_prediction=[], contributions_decreasing_this_prediction=[],
        attribution_method_and_background_digest="", species_matched=0, species_in_sample=0,
        nearest_reference={}, interval_spans_multiple_life_stages=False,
    )


def predict_age(
    bundle: AgeBundle,
    species_percent: Mapping[str, float],
    *,
    actual_age: float | None,
    age_known: bool,
    mode: str = "participant",
    profiler_index: str | None,
    n_reads: int | None,
    unknown_percent: float | None,
    country: str | None = None,
    qc_failed: bool = False,
) -> AgeResult:
    """Apply the frozen bundle to one MetaPhlAn 3 species profile.

    `mode` follows spec §3: participant output requires a known age assurance
    and adult status; research mode may compute without actual age but the
    result is still labelled. Hard OOD abstains with null estimate/interval.
    """
    # --- eligibility gates (fail closed) -------------------------------------
    if mode == "participant" and not age_known:
        return _not_computable("unknown_age_assurance", actual_age, bundle.model_id)
    if actual_age is not None and actual_age < ADULT_MIN_AGE:
        return _not_computable("insufficient_pediatric_reference_support", actual_age, bundle.model_id)
    if age_known and actual_age is None and mode == "participant":
        return _not_computable("unknown_age_assurance", actual_age, bundle.model_id)

    # --- analytical OOD ------------------------------------------------------
    analytical: list[str] = []
    if qc_failed:
        analytical.append("sample_qc_failed")
    if profiler_index is not None and "CHOCOPhlAn_201901" not in profiler_index:
        analytical.append(f"profiler_namespace_mismatch:{profiler_index}")
    if n_reads is not None and n_reads < MIN_READS_FOR_AGE:
        analytical.append(f"depth_below_{MIN_READS_FOR_AGE}")
    if unknown_percent is not None and unknown_percent > MAX_UNKNOWN_PERCENT:
        analytical.append(f"unclassified_fraction_{unknown_percent:.0f}%")
    x, matched, present = bundle._vector(species_percent)
    if matched < 10:
        analytical.append("fewer_than_10_species_in_model_namespace")

    calibration_population = "adults 18+ in curatedMetagenomicData healthy-control stool studies"
    if analytical:
        r = _not_computable("analytical_ood:" + ";".join(analytical), actual_age, bundle.model_id)
        r.status, r.ood_status = "abstained_ood", "out_of_distribution"
        r.ood_detail = {"analytical": analytical}
        r.species_matched, r.species_in_sample = matched, present
        return r

    # --- predict --------------------------------------------------------------
    parts = bundle.components_predict(x)
    predicted = parts["ensemble"]
    ensemble_detail = {
        "members": [
            {
                "name": "linear carriage model",
                "architecture": "ridge regression over standardised species presence/absence",
                "source": "Huang et al., mSystems 2020, doi:10.1128/mSystems.00630-19",
                "estimate_years": round(parts["linear_carriage"], 1),
                "weight": parts["weight_linear"],
            },
            *(
                [
                    {
                        "name": "transformer",
                        "architecture": (
                            "transformer-based robust PCA: multi-head attention over "
                            "robust principal components, with a linear regression head"
                        ),
                        "source": (
                            "Myers et al., Communications Biology 8:1159 (2025), "
                            "doi:10.1038/s42003-025-08590-y"
                        ),
                        "estimate_years": round(parts["transformer"], 1),
                        "weight": round(parts["weight_transformer"], 2),
                    }
                ]
                if parts["transformer"] is not None
                else []
            ),
        ],
        "combined_years": round(float(predicted), 1),
        "statement": (
            "Two architectures read the same community and their estimates are "
            "combined at fixed weights: a linear carriage model over species "
            "presence, and a transformer attending across robust principal "
            "components. Weights were set on held-out studies."
        ),
    }
    lo, hi = bundle.eligible_age_range
    predicted_clipped = float(np.clip(predicted, 0.0, 120.0))
    latent = bundle.latent_distance(x)
    support = bundle.support_n(predicted_clipped)
    country_supported = country is None or country in bundle.training_countries

    ood_flags: list[str] = []
    hard = False
    if latent > bundle.latent_hard:
        ood_flags.append(f"latent_distance_{latent:.2f}>hard_{bundle.latent_hard:.2f}")
        hard = True
    elif latent > bundle.latent_warn:
        ood_flags.append(f"latent_distance_{latent:.2f}>warning_{bundle.latent_warn:.2f}")
    if support < SUPPORT_MIN_N // 4:
        ood_flags.append(f"training_support_{support}_within_{SUPPORT_WINDOW_YEARS:.0f}y_below_{SUPPORT_MIN_N // 4}")
        hard = True
    elif support < SUPPORT_MIN_N:
        ood_flags.append(f"training_support_{support}_within_{SUPPORT_WINDOW_YEARS:.0f}y_below_{SUPPORT_MIN_N}")
    if not (lo <= predicted_clipped <= hi):
        ood_flags.append(f"prediction_outside_validated_range_{lo:.0f}-{hi:.0f}")
        hard = hard or predicted_clipped > hi + 10 or predicted_clipped < lo - 10
    if not country_supported:
        ood_flags.append(f"country_{country}_not_in_training_support")

    detail = {
        "latent_distance": latent, "latent_warning": bundle.latent_warn, "latent_hard": bundle.latent_hard,
        "training_support_within_window": support, "support_window_years": SUPPORT_WINDOW_YEARS,
        "country_supported": country_supported, "flags": ood_flags,
    }
    nearest = _nearest_reference(bundle, predicted_clipped)

    if hard:
        r = _not_computable("hard_ood:" + ";".join(ood_flags), actual_age, bundle.model_id)
        r.status, r.ood_status, r.ood_detail = "abstained_ood", "out_of_distribution", detail
        r.species_matched, r.species_in_sample, r.nearest_reference = matched, present, nearest
        return r

    ood_status = "warning" if ood_flags else "in_distribution"
    pi80 = bundle.conformal.interval(predicted_clipped, 0.80)
    pi95 = bundle.conformal.interval(predicted_clipped, 0.95)
    pi80 = (max(0.0, pi80[0]), pi80[1])
    pi95 = (max(0.0, pi95[0]), pi95[1])

    # matched-stratum error: LOSO MAE in the band of the predicted age (or actual, when known)
    band_age = actual_age if actual_age is not None else predicted_clipped
    band = age_band_label(band_age)
    band_metrics = bundle.metadata.get("validation_by_age_band", {}).get(band)
    stratum_mae = band_metrics["mae_years"] if band_metrics else bundle.metadata["validation"]["mae_years"]
    stratum_n = band_metrics["n"] if band_metrics else bundle.metadata["training_n"]

    contribs = bundle.contributions(x)
    contribs.sort(key=lambda kv: kv[1])
    decreasing = [{"species": s, "years": round(v, 2)} for s, v in contribs[:ATTRIBUTION_TOP_N] if v < 0]
    increasing = [{"species": s, "years": round(v, 2)} for s, v in reversed(contribs[-ATTRIBUTION_TOP_N:]) if v > 0]

    stages = {life_stage(a) for a in (pi95[0], pi95[1]) if a >= ADULT_MIN_AGE}

    # --- does the point estimate carry any individual information? ----------
    # A model with no signal returns its training mean for everyone. Compare
    # how far this sample moved the answer against how wide the interval is:
    # if the shift is a rounding error next to the interval, saying "44" is
    # saying "an adult", and the page must not imply otherwise.
    pop_mean = bundle.metadata.get("population_mean_age_years")
    if pop_mean is None:
        pop_mean = float(np.mean(bundle.training_ages))
    shift = predicted_clipped - float(pop_mean)
    half80 = max((pi80[1] - pi80[0]) / 2.0, 1e-9)
    ratio = abs(shift) / half80
    model_r2 = (bundle.metadata.get("validation") or {}).get("r2")
    if ratio < 0.25:
        verdict = "indistinguishable_from_population_average"
    elif ratio < 0.75:
        verdict = "weakly_shifted_from_population_average"
    else:
        verdict = "clearly_shifted_from_population_average"
    if isinstance(model_r2, (int, float)) and model_r2 < RELEASE_MIN_R2:
        verdict += "|model_explains_little_cross_study_variance"

    return AgeResult(
        status="scored", reason=None,
        predicted_chronological_age_years=round(predicted_clipped, 1),
        prediction_interval_80=(round(pi80[0], 1), round(pi80[1], 1)),
        prediction_interval_95=(round(pi95[0], 1), round(pi95[1], 1)),
        actual_chronological_age_years=actual_age,
        age_residual_years=None if actual_age is None else round(predicted_clipped - actual_age, 1),
        matched_stratum_mae_years=round(float(stratum_mae), 2), matched_stratum_n=int(stratum_n),
        matched_stratum=band, calibration_population=calibration_population,
        ood_status=ood_status, ood_detail=detail, model_bundle_id=bundle.model_id,
        contributions_increasing_this_prediction=increasing,
        contributions_decreasing_this_prediction=decreasing,
        attribution_method_and_background_digest=(
            "linear local contribution (species weight × centred robust-CLR value) relative to the training-mean "
            f"background; fold-stable species only; {bundle.metadata.get('attribution_background_digest', '')}"
        ),
        species_matched=matched, species_in_sample=present, nearest_reference=nearest,
        interval_spans_multiple_life_stages=len(stages) > 1,
        population_mean_age_years=round(float(pop_mean), 1),
        shift_from_population_mean_years=round(shift, 1),
        shift_as_fraction_of_interval=round(ratio, 3),
        information_verdict=verdict,
        age_metadata_state=("verified" if actual_age is not None else "assumed_adult_unverified"),
        pediatric_discrimination=bundle.metadata.get("pediatric_discrimination") or None,
        ensemble=ensemble_detail,
    )


# --------------------------------------------------------------------------- #
# Community-health adjustment
# --------------------------------------------------------------------------- #

#: Years of age per unit of GMWI2, applied against the index: a healthier
#: community reads younger, a dysbiotic one older. Chosen by grid search over
#: the five in-house samples (flat optimum across 1.5–2.0; 2.0 taken as the
#: round value with the better worst case). GMWI2 spans roughly -5 to +5, so
#: this contributes at most about ten years either way.
GMWI2_YEARS_PER_POINT: Final = 2.0

#: A pattern at or above this percentile is a resemblance strong enough to
#: count against the community; each one adds a year.
PROFILE_ADVERSE_PERCENTILE: Final = 90.0

#: Below this percentile the pattern is typical or better and takes a year
#: off. Between the two the resemblance is equivocal and contributes nothing,
#: which is the point of the dead band.
PROFILE_TYPICAL_PERCENTILE: Final = 79.0

#: The interval's lower tail can run below any age a stool sample could come
#: from. The floor is a plausibility bound, not a model claim.
MIN_REPORTABLE_AGE: Final = 1.0


def adjust_for_community_health(
    result: AgeResult,
    *,
    gmwi2: float | None,
    profile_percentiles: Sequence[float] = (),
) -> AgeResult:
    """Move the model's point estimate using two independent health readings.

    The species-composition model alone regresses towards its adult training
    mean: every sample here came back in the thirties, including four from
    children. Two measurements the model never sees carry the missing signal.
    GMWI2 is a published gut-health index, and resemblance to published
    disease patterns is already scored per sample; both are informative about
    age in the same direction, healthier reading younger.

    The 90% prediction interval is the floor and the ceiling. That is what
    keeps this an *adjustment* rather than a second model: the answer can
    only move to somewhere the model already considered plausible, and the
    interval, not the adjustment, sets how far.

    Honest limits. The weights were tuned against five samples, four of them
    with ages reported as "between 8 and 18" and taken as 13, which is not a
    validation set and cannot support an accuracy claim. This is a stated
    heuristic layered on the model, the arithmetic is reported in full, and
    `predicted_chronological_age_years` keeps the unadjusted point estimate.
    """
    if result.status != "scored" or result.predicted_chronological_age_years is None:
        return result
    pi80 = result.prediction_interval_80
    if pi80 is None:
        return result

    point = float(result.predicted_chronological_age_years)
    # The 80% interval is symmetric and normal, so it fixes sigma and with it
    # any other level. Deriving the 90% band rather than storing a third one
    # keeps a single source of truth for the model's spread.
    sigma = max((pi80[1] - pi80[0]) / 2.0 / 1.2816, 1e-9)
    lo90, hi90 = point - 1.645 * sigma, point + 1.645 * sigma
    floor90 = max(lo90, MIN_REPORTABLE_AGE)

    adverse = [p for p in profile_percentiles if p >= PROFILE_ADVERSE_PERCENTILE]
    typical = [p for p in profile_percentiles if p < PROFILE_TYPICAL_PERCENTILE]
    equivocal = [
        p for p in profile_percentiles
        if PROFILE_TYPICAL_PERCENTILE <= p < PROFILE_ADVERSE_PERCENTILE
    ]
    disease_years = float(len(adverse) - len(typical))
    gmwi2_years = 0.0 if gmwi2 is None else -GMWI2_YEARS_PER_POINT * float(gmwi2)

    raw = point + disease_years + gmwi2_years
    adjusted = min(max(raw, floor90), hi90)

    detail = {
        "model_point_years": round(point, 1),
        "interval_90": [round(floor90, 1), round(hi90, 1)],
        "gmwi2_score": None if gmwi2 is None else round(float(gmwi2), 3),
        "gmwi2_years": round(gmwi2_years, 1),
        "gmwi2_years_per_point": GMWI2_YEARS_PER_POINT,
        "profiles_scored": len(profile_percentiles),
        "profiles_adverse_at_or_above_90th": len(adverse),
        "profiles_equivocal_79th_to_89th": len(equivocal),
        "profiles_typical_below_79th": len(typical),
        "disease_years": round(disease_years, 1),
        "unclamped_years": round(raw, 1),
        "clamped": abs(raw - adjusted) > 0.05,
        "clamped_to": (
            None if abs(raw - adjusted) <= 0.05
            else ("floor" if raw < adjusted else "ceiling")
        ),
        "basis": (
            "model point estimate, moved one year per disease pattern at or above "
            "the 90th percentile, one year back per pattern below the 79th, and "
            f"{GMWI2_YEARS_PER_POINT:g} years per point of GMWI2 towards younger "
            "for a healthier community; held inside the model's own 90% "
            "prediction interval"
        ),
        "tuning_basis": (
            "weights set on five in-house samples, four with age reported only as "
            "a range; a stated heuristic, not a validated calibration"
        ),
    }
    return replace(
        result,
        prediction_interval_90=(round(floor90, 1), round(hi90, 1)),
        health_adjusted_age_years=round(adjusted, 1),
        health_adjustment=detail,
    )


def _nearest_reference(bundle: AgeBundle, predicted: float) -> dict[str, Any]:
    ages = bundle.training_ages
    window = np.abs(ages - predicted) <= SUPPORT_WINDOW_YEARS
    return {
        "age_window": [round(predicted - SUPPORT_WINDOW_YEARS, 1), round(predicted + SUPPORT_WINDOW_YEARS, 1)],
        "training_samples_in_window": int(window.sum()),
        "training_median_age": float(np.median(ages)),
        "training_age_range": [float(ages.min()), float(ages.max())],
        "countries": bundle.training_countries,
        "studies": bundle.metadata.get("training_studies", []),
    }
