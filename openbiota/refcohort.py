"""Reference cohort from curatedMetagenomicData.

v1.1.0 percentiles rested on 33 samples subsampled to 600,000 read pairs. That
is enough to rank one pathway value and not enough for what profile scoring
needs: stable percentiles at the tails, species-level taxonomy at full depth,
and references *matched* on country, age band and sex.

curatedMetagenomicData solves all three at once. It publishes uniformly
processed taxonomic profiles for ~22,600 samples across 120+ studies, with
curated subject metadata attached — country, age, sex, BMI, antibiotic use,
read depth, platform, and a study condition label. Every sample went through
one pipeline, which is the property that matters: a reference distribution
assembled from profiles computed different ways measures the pipelines, not
the biology.

Two things this module deliberately does *not* do:

* **It does not mix profiler versions.** Each cMD snapshot corresponds to one
  MetaPhlAn version. A cohort is built from a single snapshot and records it,
  because species abundances from different profiler versions are not
  interchangeable.
* **It does not silently substitute an unmatched cohort.** When a requested
  stratum is too small, it falls back to the full cohort and *says so*, so the
  confidence grade can account for it.
"""

from __future__ import annotations

import json
import math
import statistics
import time
import urllib.request
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Final

from openbiota.errors import OpenBiotaError
from openbiota.logging_util import Reporter, human_bytes, human_duration
from openbiota.net import get

#: Curated subject metadata for every cMD sample, as an R data object.
SAMPLE_METADATA_URL: Final = (
    "https://raw.githubusercontent.com/waldronlab/curatedMetagenomicData"
    "/master/data/sampleMetadata.rda"
)

#: ExperimentHub metadata, used to resolve a study's profile to a URL.
EXPERIMENTHUB_SQLITE: Final = (
    "https://experimenthub.bioconductor.org/metadata/experimenthub.sqlite3"
)

#: Default snapshot. One snapshot means one profiler version, which is the
#: whole point — see the module docstring.
DEFAULT_SNAPSHOT: Final = "2021-10-14"

#: Profiler that produced the cMD snapshots this module reads. Recorded in the
#: manifest and compared against the profiler used on the screened sample.
SNAPSHOT_PROFILER: Final = {
    "2021-03-31": "MetaPhlAn 3.0 (species-level)",
    "2021-04-02": "MetaPhlAn 3.0 (species-level)",
    "2021-10-14": "MetaPhlAn 3.0 (species-level)",
    "2022-04-13": "MetaPhlAn 3.0 (species-level)",
    "2022-10-19": "MetaPhlAn 3.0 (species-level)",
}

#: A stratum smaller than this is not used; the cohort falls back to unmatched.
MIN_STRATUM: Final = 30

#: Minimum cohort size for percentiles to be worth reporting at all.
MIN_COHORT: Final = 100

#: Age bands used for stratified matching.
AGE_BANDS: Final = ((18, 30), (30, 45), (45, 60), (60, 75), (75, 120))

#: Pseudocount for the centred log-ratio. Relative abundances contain exact
#: zeros, and log(0) is undefined; half the smallest non-zero value observed
#: is the conventional multiplicative replacement.
CLR_PSEUDO_FRACTION: Final = 0.5


def age_band(age: float | None) -> str | None:
    if age is None or not math.isfinite(age):
        return None
    for lo, hi in AGE_BANDS:
        if lo <= age < hi:
            return f"{lo}-{hi}"
    return None


# --------------------------------------------------------------------------- #
# metadata
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class CohortSample:
    """One reference sample's metadata."""

    sample_id: str
    study_name: str
    condition: str
    country: str | None
    age: float | None
    sex: str | None
    bmi: float | None
    n_reads: int | None
    platform: str | None
    antibiotics: str | None
    westernized: bool | None

    @property
    def band(self) -> str | None:
        return age_band(self.age)

    def to_json(self) -> dict[str, Any]:
        return asdict(self)


def _as_float(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _as_str(value: Any) -> str | None:
    if value is None:
        return None
    text = str(value).strip()
    return text if text and text.lower() not in ("nan", "na", "none") else None


def load_sample_metadata(cache_dir: Path, reporter: Reporter, *, refresh: bool = False) -> Any:
    """Download and parse cMD's curated subject metadata into a DataFrame."""
    try:
        import rdata
    except ImportError as exc:  # pragma: no cover — declared dependency
        raise OpenBiotaError(
            "the 'rdata' package is required to read curatedMetagenomicData; "
            "install it with `pip install rdata`"
        ) from exc

    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / "cmd_sampleMetadata.rda"
    if not path.is_file() or refresh:
        reporter.info("  downloading curatedMetagenomicData subject metadata")
        path.write_bytes(get(SAMPLE_METADATA_URL, timeout=600).body)
    parsed = rdata.read_rda(path)
    frame = next(iter(parsed.values()))
    reporter.record(f"    cMD metadata: {frame.shape[0]:,} samples x {frame.shape[1]} fields")
    return frame


def select_samples(
    frame: Any,
    *,
    condition: str = "control",
    body_site: str = "stool",
    age_categories: Sequence[str] = ("adult", "senior"),
    min_reads: int = 1_000_000,
    exclude_antibiotics: bool = True,
) -> list[CohortSample]:
    """Filter cMD metadata down to a usable cohort."""
    out: list[CohortSample] = []
    for row in frame.to_dict("records"):
        if _as_str(row.get("body_site")) != body_site:
            continue
        if condition and _as_str(row.get("study_condition")) != condition:
            continue
        category = (_as_str(row.get("age_category")) or "").lower()
        if age_categories and category not in age_categories:
            continue
        reads = _as_float(row.get("number_reads"))
        if min_reads and (reads is None or reads < min_reads):
            continue
        antibiotics = _as_str(row.get("antibiotics_current_use"))
        if exclude_antibiotics and antibiotics == "yes":
            continue
        sample_id = _as_str(row.get("sample_id"))
        study = _as_str(row.get("study_name"))
        if not sample_id or not study:
            continue
        westernized = _as_str(row.get("non_westernized"))
        out.append(
            CohortSample(
                sample_id=sample_id,
                study_name=study,
                condition=_as_str(row.get("study_condition")) or "unknown",
                country=_as_str(row.get("country")),
                age=_as_float(row.get("age")),
                sex=_as_str(row.get("gender")),
                bmi=_as_float(row.get("BMI")),
                n_reads=int(reads) if reads else None,
                platform=_as_str(row.get("sequencing_platform")),
                antibiotics=antibiotics,
                westernized=None if westernized is None else westernized == "no",
            )
        )
    return out


# --------------------------------------------------------------------------- #
# abundance profiles
# --------------------------------------------------------------------------- #


def resolve_profile_urls(
    cache_dir: Path, reporter: Reporter, *, snapshot: str = DEFAULT_SNAPSHOT
) -> dict[str, str]:
    """Map study name to its relative-abundance download URL for one snapshot."""
    import sqlite3

    cache_dir.mkdir(parents=True, exist_ok=True)
    sqlite_path = cache_dir / "experimenthub.sqlite3"
    if not sqlite_path.is_file():
        reporter.info("  downloading ExperimentHub resource index")
        sqlite_path.write_bytes(get(EXPERIMENTHUB_SQLITE, timeout=600).body)

    connection = sqlite3.connect(sqlite_path)
    query = """
        SELECT r.title, p.rdatapath, l.location_prefix
        FROM resources r
        JOIN rdatapaths p ON p.resource_id = r.id
        JOIN location_prefixes l ON l.id = r.location_prefix_id
        WHERE r.title LIKE ? AND r.title LIKE '%relative_abundance'
    """
    urls: dict[str, str] = {}
    for title, rdatapath, prefix in connection.execute(query, (f"{snapshot}.%",)):
        study = str(title).split(".")[1]
        urls[study] = f"{prefix}{rdatapath}"
    if not urls:
        raise OpenBiotaError(
            f"no relative-abundance resources found for snapshot {snapshot!r}. "
            f"Available snapshots: {', '.join(sorted(SNAPSHOT_PROFILER))}"
        )
    reporter.record(f"    snapshot {snapshot}: {len(urls)} studies with profiles")
    return urls


def download_profile(url: str, dest: Path, *, refresh: bool = False) -> Path:
    if dest.is_file() and not refresh:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "openbiota"})
    with urllib.request.urlopen(request, timeout=600) as response:
        dest.write_bytes(response.read())
    return dest


def parse_profile(path: Path) -> tuple[list[str], list[str], list[list[float]]]:
    """Read a cMD relative-abundance matrix.

    Returns ``(taxa, samples, values)`` restricted to species-level rows.
    Species is the only level a profile may use — see PANELS/PROFILES docs and
    the conflicted-genus discipline in the spec.
    """
    import rdata

    parsed = rdata.read_rda(path)
    array = next(iter(parsed.values()))
    taxa_all = [str(x) for x in array.coords["dim_0"].values]
    samples = [str(x) for x in array.coords["dim_1"].values]
    values = array.values

    keep = [
        index
        for index, name in enumerate(taxa_all)
        if "|s__" in name and "|t__" not in name
    ]
    taxa = [species_name(taxa_all[i]) for i in keep]
    matrix = [[float(values[i][j]) for j in range(len(samples))] for i in keep]
    return taxa, samples, matrix


def species_name(lineage: str) -> str:
    """Extract the bare species name from a MetaPhlAn lineage string."""
    for part in reversed(lineage.split("|")):
        if part.startswith("s__"):
            return part[3:]
    return lineage


# --------------------------------------------------------------------------- #
# the assembled cohort
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class ReferenceCohort:
    """Species abundances plus metadata for a set of reference samples."""

    snapshot: str
    profiler: str
    taxa: list[str]
    sample_ids: list[str]
    #: taxon index -> sample index -> relative abundance (percent)
    abundance: list[list[float]]
    metadata: dict[str, CohortSample]
    manifest_id: str = ""
    built_at: str = ""
    _frame_cache: list[int] | None = field(default=None, repr=False, compare=False)
    _frame_mean_log_cache: list[float] | None = field(default=None, repr=False, compare=False)
    _index_cache: dict[str, int] | None = field(default=None, repr=False, compare=False)
    _taxon_refs_cache: dict[str, TaxonReference] | None = field(
        default=None, repr=False, compare=False
    )

    @property
    def n_samples(self) -> int:
        return len(self.sample_ids)

    @property
    def n_taxa(self) -> int:
        return len(self.taxa)

    @property
    def frame(self) -> list[int]:
        """Prevalent-core denominator for the CLR; see `CLR_FRAME_MIN_PREVALENCE`."""
        cached = getattr(self, "_frame_cache", None)
        if cached is None:
            cached = prevalent_frame(self.abundance)
            object.__setattr__(self, "_frame_cache", cached)
        return cached

    @property
    def frame_mean_log(self) -> list[float]:
        """Per-sample geometric-mean denominator (mean log over the frame).

        This is the same denominator ``clr`` uses; caching it lets group
        aggregates (genus sums, carrier sums) be placed on the CLR scale
        without re-walking every taxon of every sample for every group.
        """
        cached = getattr(self, "_frame_mean_log_cache", None)
        if cached is None:
            frame = self.frame or list(range(self.n_taxa))
            cached = []
            for j in range(self.n_samples):
                logs = [
                    math.log(v if v > 0 else CLR_FLOOR_PERCENT)
                    for v in (self.abundance[i][j] for i in frame)
                ]
                cached.append(sum(logs) / len(logs) if logs else 0.0)
            object.__setattr__(self, "_frame_mean_log_cache", cached)
        return cached

    def index_of(self, taxon: str) -> int | None:
        cached = getattr(self, "_index_cache", None)
        if cached is None or len(cached) != len(self.taxa):
            cached = {name: i for i, name in enumerate(self.taxa)}
            object.__setattr__(self, "_index_cache", cached)
        return cached.get(taxon)

    def composition(self) -> dict[str, Any]:
        """Cohort make-up, for the manifest and the report."""
        from collections import Counter

        samples = [self.metadata[s] for s in self.sample_ids if s in self.metadata]
        return {
            "n_samples": len(samples),
            "n_studies": len({s.study_name for s in samples}),
            "countries": dict(Counter(s.country or "unknown" for s in samples).most_common()),
            "age_bands": dict(Counter(s.band or "unknown" for s in samples).most_common()),
            "sex": dict(Counter(s.sex or "unknown" for s in samples).most_common()),
            "conditions": dict(Counter(s.condition for s in samples).most_common()),
            "median_reads": (
                statistics.median([s.n_reads for s in samples if s.n_reads])
                if any(s.n_reads for s in samples)
                else None
            ),
        }

    def subset(self, sample_ids: Iterable[str]) -> ReferenceCohort:
        wanted = [s for s in sample_ids if s in set(self.sample_ids)]
        positions = [self.sample_ids.index(s) for s in wanted]
        return ReferenceCohort(
            snapshot=self.snapshot,
            profiler=self.profiler,
            taxa=list(self.taxa),
            sample_ids=wanted,
            abundance=[[row[p] for p in positions] for row in self.abundance],
            metadata=self.metadata,
            manifest_id=self.manifest_id,
            built_at=self.built_at,
        )

    def match(
        self,
        *,
        country: str | None = None,
        band: str | None = None,
        sex: str | None = None,
        min_stratum: int = MIN_STRATUM,
    ) -> tuple[ReferenceCohort, dict[str, Any]]:
        """Stratified subset, falling back to the full cohort when too small.

        The fallback is reported rather than hidden — an unmatched reference is
        a legitimate answer, but it has to feed into the confidence grade.
        """
        criteria = {"country": country, "age_band": band, "sex": sex}
        active = {k: v for k, v in criteria.items() if v}
        if not active:
            return self, {"matched": False, "reason": "no matching criteria requested",
                          "criteria": {}, "n": self.n_samples}

        def keep(sample: CohortSample) -> bool:
            if country and sample.country != country:
                return False
            if band and sample.band != band:
                return False
            return not (sex and sample.sex != sex)

        chosen = [s for s in self.sample_ids if s in self.metadata and keep(self.metadata[s])]
        if len(chosen) < min_stratum:
            return self, {
                "matched": False,
                "reason": (
                    f"stratum {active} holds {len(chosen)} samples, below the minimum of "
                    f"{min_stratum}; using the full cohort instead"
                ),
                "criteria": active,
                "n": self.n_samples,
                "stratum_n": len(chosen),
            }
        return self.subset(chosen), {
            "matched": True,
            "reason": f"matched on {', '.join(active)}",
            "criteria": active,
            "n": len(chosen),
        }

    # -- persistence ------------------------------------------------------- #

    def manifest(self) -> dict[str, Any]:
        return {
            "manifest_id": self.manifest_id,
            "snapshot": self.snapshot,
            "profiler": self.profiler,
            "built_at": self.built_at,
            "n_samples": self.n_samples,
            "n_taxa": self.n_taxa,
            "composition": self.composition(),
            "source": "curatedMetagenomicData",
            "citation": (
                "Pasolli E. et al., Nat Methods 14:1023-1024 (2017), doi:10.1038/nmeth.4468"
            ),
        }

    def save(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "manifest": self.manifest(),
            "taxa": self.taxa,
            "sample_ids": self.sample_ids,
            "abundance": self.abundance,
            "metadata": {k: v.to_json() for k, v in self.metadata.items()
                         if k in set(self.sample_ids)},
        }
        path.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> ReferenceCohort:
        """The cohort at `path`, from a fingerprinted cache when the file is unchanged.

        The JSON is 14 MB, almost all of it the abundance matrix, and parsing
        it costs ~0.4 s on every report run. The pickle carries the same
        object and loads in a fraction of that.
        """
        from openbiota.fastcache import cached_load

        return cached_load("taxonomic_cohort", [path], lambda: cls._load_uncached(path))

    @classmethod
    def _load_uncached(cls, path: Path) -> ReferenceCohort:
        payload = json.loads(path.read_text(encoding="utf-8"))
        manifest = payload["manifest"]
        return cls(
            snapshot=manifest["snapshot"],
            profiler=manifest["profiler"],
            taxa=payload["taxa"],
            sample_ids=payload["sample_ids"],
            abundance=payload["abundance"],
            metadata={k: CohortSample(**v) for k, v in payload["metadata"].items()},
            manifest_id=manifest.get("manifest_id", ""),
            built_at=manifest.get("built_at", ""),
        )


def build_cohort(
    *,
    cache_dir: Path,
    reporter: Reporter,
    snapshot: str = DEFAULT_SNAPSHOT,
    condition: str = "control",
    max_samples: int | None = None,
    min_reads: int = 1_000_000,
    refresh: bool = False,
) -> ReferenceCohort:
    """Download and assemble a reference cohort from one cMD snapshot."""
    started = time.monotonic()
    frame = load_sample_metadata(cache_dir, reporter, refresh=refresh)
    wanted = select_samples(frame, condition=condition, min_reads=min_reads)
    reporter.info(
        f"  {len(wanted):,} samples match (condition={condition}, stool, adult, "
        f">={min_reads:,} reads, no current antibiotics)"
    )
    if not wanted:
        raise OpenBiotaError(f"no cMD samples matched condition={condition!r}")

    by_study: dict[str, list[CohortSample]] = {}
    for sample in wanted:
        by_study.setdefault(sample.study_name, []).append(sample)

    urls = resolve_profile_urls(cache_dir, reporter, snapshot=snapshot)
    available = [s for s in by_study if s in urls]
    missing = sorted(set(by_study) - set(urls))
    if missing:
        reporter.record(
            f"    {len(missing)} studies have no profile in this snapshot: "
            f"{', '.join(missing[:6])}{' …' if len(missing) > 6 else ''}"
        )

    taxa_union: dict[str, int] = {}
    columns: list[tuple[str, dict[str, float]]] = []
    downloaded = 0

    for index, study in enumerate(sorted(available), start=1):
        path = cache_dir / "profiles" / snapshot / f"{study}.rda"
        try:
            download_profile(urls[study], path, refresh=refresh)
        except Exception as exc:  # noqa: BLE001 — one bad study must not stop the build
            reporter.warn(f"  {study}: profile download failed ({exc}); skipped")
            continue
        downloaded += path.stat().st_size
        try:
            taxa, samples, matrix = parse_profile(path)
        except Exception as exc:  # noqa: BLE001
            reporter.warn(f"  {study}: profile unreadable ({exc}); skipped")
            continue

        wanted_ids = {s.sample_id for s in by_study[study]}
        for position, sample_id in enumerate(samples):
            if sample_id not in wanted_ids:
                continue
            column = {
                taxa[i]: matrix[i][position]
                for i in range(len(taxa))
                if matrix[i][position] > 0
            }
            for taxon in column:
                taxa_union.setdefault(taxon, len(taxa_union))
            columns.append((sample_id, column))

        if index % 15 == 0 or index == len(available):
            reporter.info(
                f"  cohort {index}/{len(available)} studies, {len(columns):,} samples, "
                f"{human_bytes(downloaded)} downloaded, "
                f"{human_duration(time.monotonic() - started)}"
            )
        if max_samples and len(columns) >= max_samples:
            reporter.info(f"  reached --max-samples {max_samples:,}; stopping")
            break

    if len(columns) < MIN_COHORT:
        raise OpenBiotaError(
            f"assembled only {len(columns)} samples; at least {MIN_COHORT} are needed for "
            "percentiles to mean anything"
        )

    taxa = sorted(taxa_union)
    positions = {taxon: i for i, taxon in enumerate(taxa)}
    abundance = [[0.0] * len(columns) for _ in taxa]
    for column_index, (_, column) in enumerate(columns):
        for taxon, value in column.items():
            abundance[positions[taxon]][column_index] = value

    metadata = {s.sample_id: s for s in wanted}
    manifest_id = _manifest_id(snapshot, [sample_id for sample_id, _ in columns])
    cohort = ReferenceCohort(
        snapshot=snapshot,
        profiler=SNAPSHOT_PROFILER.get(snapshot, "unknown"),
        taxa=taxa,
        sample_ids=[sample_id for sample_id, _ in columns],
        abundance=abundance,
        metadata={k: v for k, v in metadata.items() if k in {c[0] for c in columns}},
        manifest_id=manifest_id,
        built_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
    )
    reporter.ok(
        f"reference cohort: {cohort.n_samples:,} samples, {cohort.n_taxa:,} species, "
        f"manifest {manifest_id}"
    )
    return cohort


def _manifest_id(snapshot: str, sample_ids: Sequence[str]) -> str:
    import hashlib

    payload = snapshot + "|" + "|".join(sorted(sample_ids))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]


# --------------------------------------------------------------------------- #
# compositional transform
# --------------------------------------------------------------------------- #


#: Fixed zero-replacement floor, as a relative abundance percentage.
#:
#: This MUST be a global constant rather than derived per sample. Deriving it
#: from each sample's own smallest detected value — the textbook
#: "half the minimum" rule — makes the floor depth-dependent: a deeply
#: sequenced sample detects rarer taxa, gets a lower floor, and therefore maps
#: its zeros to more negative CLR values than a shallow one. For a taxon that
#: is absent in most samples, that difference *is* the signal, so the score
#: ends up tracking sequencing depth. Measured on the CRC validation set,
#: switching to a fixed floor is the difference between an uninterpretable
#: score and a working one.
#:
#: 1e-4 percent sits just below MetaPhlAn's practical detection limit at
#: typical depth, so a real detection always lands above it.
CLR_FLOOR_PERCENT: Final = 1e-4


#: Taxa at or above this prevalence in the reference form the CLR denominator.
#:
#: The textbook CLR divides by the geometric mean of *every* component. For a
#: taxon near the detection floor, that makes its transformed value depend on
#: how many rare taxa happened to be detected in that sample — noise that
#: has nothing to do with the taxon itself. Dividing instead by the geometric
#: mean of a stable core of prevalent taxa (a "reference frame") removes it.
#: Measured on the CRC validation set this recovers a third of the gap between
#: whole-composition CLR and the theoretical ceiling for the features.
CLR_FRAME_MIN_PREVALENCE: Final = 0.30


def clr(
    values: Sequence[float],
    *,
    pseudocount: float = CLR_FLOOR_PERCENT,
    frame: Sequence[int] | None = None,
) -> list[float]:
    """Centred log-ratio transform with a fixed floor and optional reference frame.

    Relative abundances are compositional: they sum to a constant, so an
    apparent increase in one taxon can be produced entirely by decreases
    elsewhere, and raw sums or differences of them are not valid. The CLR maps
    the composition into a real space where they are.

    ``pseudocount`` replaces zeros and is a fixed constant for the reason given
    at `CLR_FLOOR_PERCENT`. ``frame`` optionally restricts the geometric-mean
    denominator to a stable set of component indices; see
    `CLR_FRAME_MIN_PREVALENCE`. Sample and reference must use the same floor
    and the same frame or the two sides of the comparison are on different
    scales.
    """
    if not values:
        return []
    logs = [math.log(v if v > 0 else pseudocount) for v in values]
    denominator = [logs[i] for i in frame] if frame else logs
    if not denominator:
        denominator = logs
    mean_log = statistics.fmean(denominator)
    return [lp - mean_log for lp in logs]


def prevalent_frame(
    abundance: Sequence[Sequence[float]], *, min_prevalence: float = CLR_FRAME_MIN_PREVALENCE
) -> list[int]:
    """Indices of taxa present in at least ``min_prevalence`` of samples."""
    if not abundance or not abundance[0]:
        return []
    n = len(abundance[0])
    frame = [
        i
        for i, row in enumerate(abundance)
        if sum(1 for v in row if v > 0) / n >= min_prevalence
    ]
    # A frame needs enough members to be stable; below that, fall back to the
    # whole composition rather than a handful of taxa.
    return frame if len(frame) >= 20 else []


@dataclass(frozen=True, slots=True)
class TaxonReference:
    """Robust location and scale for one taxon across the reference cohort."""

    taxon: str
    median: float
    mad: float
    n: int
    prevalence: float
    #: CLR values, retained so a percentile can be taken without rescaling
    values: tuple[float, ...] = ()

    def robust_z(self, value: float) -> float | None:
        """MAD-based z with the same degenerate-spread fallback as scoring.

        Delegates so that the cohort summary and the scoring engine cannot
        drift apart on the one piece of arithmetic they share.
        """
        from openbiota.scoring import FeatureReference

        return FeatureReference(
            feature=self.taxon,
            median=self.median,
            mad=self.mad,
            n=self.n,
            source="cohort",
            values=self.values,
        ).robust_z(value)

    def percentile_of(self, value: float) -> float | None:
        if not self.values:
            return None
        below = sum(1 for v in self.values if v < value)
        return below / len(self.values) * 100.0

    def to_json(self) -> dict[str, Any]:
        return {
            "taxon": self.taxon,
            "clr_median": round(self.median, 4),
            "clr_mad": round(self.mad, 4),
            "n": self.n,
            "prevalence": round(self.prevalence, 4),
        }


def build_taxon_references(cohort: ReferenceCohort) -> dict[str, TaxonReference]:
    """CLR-transform the cohort and summarise each taxon robustly.

    The CLR is computed per sample across all taxa — that is what makes it a
    compositional transform — and the summary is then taken per taxon across
    samples. The denominator uses the cohort's prevalent-core frame so that a
    sample scored later against these references is transformed identically.
    """
    cached = getattr(cohort, "_taxon_refs_cache", None)
    if cached is not None:
        return cached
    n_samples = cohort.n_samples
    frame = cohort.frame
    clr_columns: list[list[float]] = []
    for j in range(n_samples):
        column = [cohort.abundance[i][j] for i in range(cohort.n_taxa)]
        clr_columns.append(clr(column, frame=frame))

    references: dict[str, TaxonReference] = {}
    for i, taxon in enumerate(cohort.taxa):
        series = [clr_columns[j][i] for j in range(n_samples)]
        median = statistics.median(series)
        mad = statistics.median([abs(v - median) for v in series])
        prevalence = sum(1 for j in range(n_samples) if cohort.abundance[i][j] > 0) / max(
            n_samples, 1
        )
        references[taxon] = TaxonReference(
            taxon=taxon,
            median=median,
            mad=mad,
            n=n_samples,
            prevalence=prevalence,
            values=tuple(series),
        )
    # The cohort object is immutable once built, so its references are too;
    # the similarity stage and the community overview both ask for them.
    object.__setattr__(cohort, "_taxon_refs_cache", references)
    return references


def sample_clr(
    abundances: Mapping[str, float], cohort: ReferenceCohort
) -> dict[str, float]:
    """CLR-transform one sample onto the cohort's taxon axis and frame.

    The sample is projected onto the reference taxon list so that both sides
    of the comparison span the same components, and uses the same denominator
    frame — a CLR computed over a different set of parts is not comparable.
    """
    vector = [max(0.0, abundances.get(taxon, 0.0)) for taxon in cohort.taxa]
    transformed = clr(vector, frame=cohort.frame)
    return dict(zip(cohort.taxa, transformed, strict=True))


# --------------------------------------------------------------------------- #
# stability curve (spec 7.5)
# --------------------------------------------------------------------------- #


def stability_curve(
    references: dict[str, TaxonReference],
    *,
    taxa: Sequence[str],
    sizes: Sequence[int] = (25, 50, 100, 200, 400, 800),
    replicates: int = 40,
    seed: int = 20260904,
) -> list[dict[str, Any]]:
    """How much does a percentile estimate move as the cohort grows?

    Bootstraps sub-cohorts of increasing size and measures the spread of the
    median estimate. Justifies a cohort-size target empirically instead of by
    assertion.
    """
    import random

    rng = random.Random(seed)
    usable = [t for t in taxa if t in references and references[t].values]
    if not usable:
        return []

    out: list[dict[str, Any]] = []
    full_n = len(references[usable[0]].values)
    for size in sizes:
        if size > full_n:
            continue
        spreads: list[float] = []
        for taxon in usable:
            series = references[taxon].values
            medians = [
                statistics.median(rng.sample(series, size)) for _ in range(replicates)
            ]
            if len(medians) > 1:
                spreads.append(statistics.stdev(medians))
        if spreads:
            out.append(
                {
                    "cohort_size": size,
                    "median_estimate_sd": round(statistics.fmean(spreads), 4),
                    "taxa_assessed": len(usable),
                    "replicates": replicates,
                }
            )
    return out
