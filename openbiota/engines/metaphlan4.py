"""Extended catalogue lane: MetaPhlAn 4 with species-level genome bins.

Two profilers, two jobs
-----------------------
The scoring lane stays on MetaPhlAn 3 (``openbiota.engines.metaphlan``) because
every reference it is compared against — the curatedMetagenomicData cohort,
the GMWI2 anchor, the published disease signatures — was built on MetaPhlAn 3
species. A percentile is only meaningful when sample and reference were named
the same way.

This lane exists to answer a different question: *what is in the sample*, as
completely as current references allow. MetaPhlAn 4's CHOCOPhlAnSGB database
(Blanco-Míguez et al., Nat Biotechnol 2023) names ~26,000 species-level genome
bins (SGBs) built from 1.01 million genomes and metagenome-assembled genomes,
against ~13,500 species in the MetaPhlAn 3 catalogue. Around 4,900 of the SGBs
are "unknown" (uSGB): organisms with no cultured representative and no formal
name, which MetaPhlAn 3 could not report at all. It also carries the 2020
Lactobacillus split, Phocaeicola, Segatella and the other renamings.

The output of this lane is reported, never scored. It fills the species
inventory, the archaea/eukaryote counts, and the "organisms without a name"
figure; the disease profiles do not read it.

Same host-filtered reads
------------------------
The lane reuses the non-host FASTQ produced for the scoring lane, so the two
profilers see identical input and the host fraction is measured once.
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

from openbiota.engines.metaphlan import _env_for, _fingerprint, _metaphlan_version
from openbiota.errors import DependencyError, OpenBiotaError
from openbiota.logging_util import Reporter, human_duration

#: Pinned SGB database. Changing it invalidates every cached extended profile.
PINNED_INDEX4: Final = "mpa_vJun23_CHOCOPhlAnSGB_202403"
CACHE_VERSION: Final = 1

#: MetaPhlAn 3 species name -> current taxonomic name, for the well-known
#: renames that would otherwise make the two inventories look more different
#: than they are. The Jun23 SGB catalogue applies some of these (Phocaeicola,
#: Segatella, the Lactobacillus split) and keeps the old name for others
#: (Eubacterium rectale, Ruminococcus torques), so a name is treated as
#: shared when *either* form is present. Not exhaustive: everything else is
#: matched exactly or reported as "named only in one".
RENAMED: Final[dict[str, str]] = {
    "Bacteroides_vulgatus": "Phocaeicola_vulgatus",
    "Bacteroides_dorei": "Phocaeicola_dorei",
    "Bacteroides_plebeius": "Phocaeicola_plebeius",
    "Bacteroides_coprocola": "Phocaeicola_coprocola",
    "Bacteroides_massiliensis": "Phocaeicola_massiliensis",
    "Bacteroides_coprophilus": "Phocaeicola_coprophilus",
    "Bacteroides_sartorii": "Phocaeicola_sartorii",
    "Prevotella_copri": "Segatella_copri",
    "Eubacterium_rectale": "Agathobacter_rectalis",
    "Eubacterium_hallii": "Anaerobutyricum_hallii",
    "Eubacterium_eligens": "Lachnospira_eligens",
    "Ruminococcus_gnavus": "Mediterraneibacter_gnavus",
    "Ruminococcus_torques": "Mediterraneibacter_torques",
    "Ruminococcus_obeum": "Blautia_obeum",
    "Ruminococcus_lactaris": "Mediterraneibacter_lactaris",
    "Clostridium_bolteae": "Enterocloster_bolteae",
    "Clostridium_clostridioforme": "Enterocloster_clostridioformis",
    "Clostridium_asparagiforme": "Enterocloster_asparagiformis",
    "Clostridium_citroniae": "Enterocloster_citroniae",
    "Clostridium_hathewayi": "Hungatella_hathewayi",
    "Clostridium_symbiosum": "Clostridium_symbiosum",
    "Clostridium_leptum": "Clostridium_leptum",
    "Lactobacillus_rhamnosus": "Lacticaseibacillus_rhamnosus",
    "Lactobacillus_casei": "Lacticaseibacillus_casei",
    "Lactobacillus_paracasei": "Lacticaseibacillus_paracasei",
    # Bacillus coagulans moved to Weizmannia in 2020 (Gupta et al.), and the
    # newer catalogue uses the current name. Verified against this sample:
    # the older lane reports Bacillus_coagulans where the newer reports
    # Weizmannia_coagulans, and never both.
    "Bacillus_coagulans": "Weizmannia_coagulans",
    "Lactobacillus_plantarum": "Lactiplantibacillus_plantarum",
    "Lactobacillus_reuteri": "Limosilactobacillus_reuteri",
    "Lactobacillus_fermentum": "Limosilactobacillus_fermentum",
    "Lactobacillus_mucosae": "Limosilactobacillus_mucosae",
    "Lactobacillus_salivarius": "Ligilactobacillus_salivarius",
    "Lactobacillus_ruminis": "Ligilactobacillus_ruminis",
    "Lactobacillus_brevis": "Levilactobacillus_brevis",
    "Lactobacillus_sakei": "Latilactobacillus_sakei",
    "Lactobacillus_curvatus": "Latilactobacillus_curvatus",
    "Eubacterium_siraeum": "Eubacterium_siraeum",
    "Eubacterium_ventriosum": "Eubacterium_ventriosum",
    "Coprococcus_comes": "Coprococcus_comes",
    "Dorea_formicigenerans": "Dorea_formicigenerans",
    "Roseburia_faecis": "Roseburia_faecis",
    "Streptococcus_thermophilus": "Streptococcus_thermophilus",
    "Bifidobacterium_adolescentis": "Bifidobacterium_adolescentis",
}


@dataclass(frozen=True, slots=True)
class Metaphlan4Tools:
    metaphlan: str
    version: str


def locate_tools4(*, explicit: str | None = None) -> Metaphlan4Tools:
    """Find a MetaPhlAn 4 executable: explicit path, ``$OPENBIOTA_METAPHLAN4``,
    the project's ``.venv-mpa4``, then ``PATH``."""
    candidates: list[str] = []
    if explicit:
        candidates.append(explicit)
    if os.environ.get("OPENBIOTA_METAPHLAN4"):
        candidates.append(os.environ["OPENBIOTA_METAPHLAN4"])
    repo_root = Path(__file__).resolve().parents[2]
    candidates.append(str(repo_root / ".venv-mpa4" / "bin" / "metaphlan"))
    on_path = shutil.which("metaphlan")
    if on_path:
        candidates.append(on_path)
    for candidate in candidates:
        if not Path(candidate).is_file():
            continue
        version = _metaphlan_version(candidate)
        if version.startswith("4."):
            return Metaphlan4Tools(metaphlan=candidate, version=version)
    raise DependencyError(
        "MetaPhlAn 4 not found. Create the dedicated environment with\n"
        "    python3 -m venv .venv-mpa4 && .venv-mpa4/bin/pip install 'metaphlan>=4.1'\n"
        "or set OPENBIOTA_METAPHLAN4 to a metaphlan 4 executable. The extended catalogue "
        "is optional; the report runs without it."
    )


def database4_ready(db_dir: Path, index: str = PINNED_INDEX4) -> bool:
    bt2 = [db_dir / f"{index}.{n}.bt2l" for n in ("1", "2", "3", "4", "rev.1", "rev.2")]
    small = [db_dir / f"{index}.{n}.bt2" for n in ("1", "2", "3", "4", "rev.1", "rev.2")]
    return (all(p.is_file() for p in bt2) or all(p.is_file() for p in small)) and (
        db_dir / f"{index}.pkl"
    ).is_file()


# --------------------------------------------------------------------------- #
# result
# --------------------------------------------------------------------------- #


def names_for_bridge(mpa3_names: set[str]) -> set[str]:
    """Every name under which a MetaPhlAn 3 species might appear in MetaPhlAn 4."""
    out: set[str] = set()
    for name in mpa3_names:
        out.add(name)
        out.add(RENAMED.get(name, name))
    return out


@dataclass(frozen=True, slots=True)
class SGBRow:
    """One SGB-level entry from the extended catalogue."""

    sgb: str
    species: str
    genus: str
    kingdom: str
    percent: float
    #: True for uSGBs: no cultured isolate, no formal name (``GGB…_SGB…``).
    unnamed: bool

    def to_json(self) -> dict[str, Any]:
        return {
            "sgb": self.sgb,
            "species": self.species,
            "genus": self.genus,
            "kingdom": self.kingdom,
            "percent": round(self.percent, 4),
            "unnamed": self.unnamed,
        }


@dataclass(frozen=True, slots=True)
class ExtendedCatalogue:
    """One sample's MetaPhlAn 4 SGB profile, parsed and summarised."""

    index: str
    profiler_version: str
    clades: dict[str, float]
    unclassified_percent: float
    n_reads_processed: int
    elapsed_s: float
    cached: bool
    command: tuple[str, ...] = field(default=())

    def level(self, prefix: str) -> dict[str, float]:
        out: dict[str, float] = {}
        for lineage, value in self.clades.items():
            last = lineage.rsplit("|", 1)[-1]
            if last.startswith(prefix) and self._is_terminal_for(lineage, prefix):
                out[last[len(prefix):]] = value
        return out

    @staticmethod
    def _is_terminal_for(lineage: str, prefix: str) -> bool:
        # A species line is "...|s__X" (no t__); an SGB line ends in "|t__SGB123".
        last = lineage.rsplit("|", 1)[-1]
        return last.startswith(prefix)

    @property
    def species(self) -> dict[str, float]:
        return self.level("s__")

    @property
    def genera(self) -> dict[str, float]:
        return self.level("g__")

    @property
    def phyla(self) -> dict[str, float]:
        return self.level("p__")

    @property
    def kingdoms(self) -> dict[str, float]:
        return self.level("k__")

    def rows(self) -> list[SGBRow]:
        """SGB-level rows with their species, genus and kingdom, largest first."""
        out: list[SGBRow] = []
        for lineage, value in self.clades.items():
            if "|t__" not in lineage or value <= 0:
                continue
            parts = dict(p.split("__", 1) for p in lineage.split("|") if "__" in p)
            sgb = parts.get("t", "")
            species = parts.get("s", "")
            out.append(
                SGBRow(
                    sgb=sgb,
                    species=species,
                    genus=parts.get("g", ""),
                    kingdom=parts.get("k", ""),
                    percent=value,
                    unnamed=species.startswith("GGB") or "_SGB" in species,
                )
            )
        return sorted(out, key=lambda r: -r.percent)

    @property
    def n_sgbs(self) -> int:
        return sum(1 for r in self.rows())

    @property
    def n_species(self) -> int:
        return sum(1 for v in self.species.values() if v > 0)

    @property
    def n_genera(self) -> int:
        return sum(1 for v in self.genera.values() if v > 0)

    @property
    def n_unnamed(self) -> int:
        return sum(1 for r in self.rows() if r.unnamed)

    @property
    def unnamed_percent(self) -> float:
        return sum(r.percent for r in self.rows() if r.unnamed)

    def kingdom_counts(self) -> dict[str, int]:
        counts: dict[str, int] = {}
        for r in self.rows():
            counts[r.kingdom] = counts.get(r.kingdom, 0) + 1
        return counts

    def bridge(self, mpa3_species: dict[str, float]) -> dict[str, Any]:
        """Compare with the scoring lane's species list.

        Returns counts of species named by both, by MetaPhlAn 4 only (new
        catalogue entries and uSGBs) and by MetaPhlAn 3 only (mostly names
        that were split or merged).
        """
        mpa4 = {k for k, v in self.species.items() if v > 0}
        mpa3 = {k for k, v in mpa3_species.items() if v > 0}
        translated = names_for_bridge(mpa3)
        both = mpa4 & translated
        only4 = mpa4 - translated
        only3 = {k for k in mpa3 if not ({k, RENAMED.get(k, k)} & mpa4)}
        return {
            "named_by_both": len(both),
            "named_only_by_metaphlan4": len(only4),
            "named_only_by_metaphlan3": len(only3),
            "metaphlan4_only_examples": sorted(only4, key=lambda k: -self.species[k])[:12],
            "metaphlan3_only_examples": sorted(only3, key=lambda k: -mpa3_species[k])[:12],
        }

    def to_json(self, *, mpa3_species: dict[str, float] | None = None) -> dict[str, Any]:
        rows = self.rows()
        return {
            "engine": "metaphlan4",
            "profiler_version": self.profiler_version,
            "database": self.index,
            "unclassified_percent": round(self.unclassified_percent, 3),
            "n_reads_processed": self.n_reads_processed,
            "n_sgbs_detected": len(rows),
            "n_species_detected": self.n_species,
            "n_genera_detected": self.n_genera,
            "n_unnamed_sgbs": self.n_unnamed,
            "unnamed_sgb_percent": round(self.unnamed_percent, 3),
            "kingdom_counts": self.kingdom_counts(),
            "kingdom_percent": {k: round(v, 3) for k, v in self.kingdoms.items()},
            "phyla": {k: round(v, 3) for k, v in sorted(self.phyla.items(), key=lambda kv: -kv[1])},
            "sgbs": [r.to_json() for r in rows],
            "bridge_to_scoring_lane": None if mpa3_species is None else self.bridge(mpa3_species),
            "cached": self.cached,
            "elapsed_s": round(self.elapsed_s, 1),
            "what_this_is": (
                "The most complete inventory current references allow: MetaPhlAn 4 against "
                "CHOCOPhlAnSGB, which names ~26,000 species-level genome bins including "
                "~4,900 with no cultured representative. Reported for completeness; the "
                "disease profiles and percentiles use the MetaPhlAn 3 lane so that sample "
                "and reference are named the same way."
            ),
        }


# --------------------------------------------------------------------------- #
# running
# --------------------------------------------------------------------------- #


def parse_profile4(path: Path) -> tuple[dict[str, float], float, int]:
    """Parse a MetaPhlAn 4 profile into (clades, unclassified_percent, n_reads)."""
    clades: dict[str, float] = {}
    unclassified = 0.0
    n_reads = 0
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            if "reads processed" in line:
                digits = "".join(ch for ch in line if ch.isdigit())
                n_reads = int(digits) if digits else 0
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        name = parts[0]
        try:
            value = float(parts[2])
        except ValueError:
            continue
        if name in ("UNCLASSIFIED", "UNKNOWN"):
            unclassified = value
            continue
        clades[name] = value
    return clades, unclassified, n_reads


def run_metaphlan4(
    tools: Metaphlan4Tools,
    *,
    r1: Path,
    r2: Path | None,
    db_dir: Path,
    work_dir: Path,
    reporter: Reporter,
    threads: int = 8,
    index: str = PINNED_INDEX4,
    sam_out: Path | None = None,
) -> ExtendedCatalogue:
    """Profile with MetaPhlAn 4; cache the Bowtie2 alignment and the profile.

    `sam_out` asks for the marker alignment as SAM as well, which is what
    strain typing (`sample2markers`) consumes. MetaPhlAn only writes SAM
    while it is running Bowtie2, so when the SAM is wanted and missing the
    alignment is redone from the reads even if a cached `bowtie2out`
    exists: the cache is what let the strain lane silently not run.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    bowtie2out = work_dir / f"metaphlan4.{index}.bowtie2.bz2"
    profile_path = work_dir / f"metaphlan4.{index}.tsv"
    stamp = work_dir / f"metaphlan4.{index}.json"
    need_sam = sam_out is not None and not (sam_out.is_file() and sam_out.stat().st_size > 0)

    fingerprint = _fingerprint(r1, r2)
    if stamp.is_file() and profile_path.is_file() and not need_sam:
        meta = json.loads(stamp.read_text())
        if meta.get("cache_version") == CACHE_VERSION and meta.get("inputs") == fingerprint:
            clades, unclassified, n_reads = parse_profile4(profile_path)
            reporter.record(f"    extended catalogue: cached ({profile_path.name})")
            return ExtendedCatalogue(
                index=index, profiler_version=meta.get("version", tools.version),
                clades=clades, unclassified_percent=unclassified, n_reads_processed=n_reads,
                elapsed_s=meta.get("elapsed_s", 0.0), cached=True,
                command=tuple(meta.get("command", ())),
            )

    if not database4_ready(db_dir, index):
        raise DependencyError(
            f"MetaPhlAn 4 database {index} not found in {db_dir}. Download the tarball from "
            "http://cmprod1.cibio.unitn.it/biobakery4/metaphlan_databases/ and unpack it "
            "there, or run `openbiota install-mpa4`."
        )

    inputs = str(r1) if r2 is None else f"{r1},{r2}"
    command = [
        tools.metaphlan, inputs,
        "--input_type", "fastq",
        "--bowtie2db", str(db_dir),
        "--index", index,
        "--nproc", str(threads),
        "--unclassified_estimation",
        "-o", str(profile_path),
    ]
    if bowtie2out.is_file() and not need_sam:
        command = [
            tools.metaphlan, str(bowtie2out), "--input_type", "bowtie2out",
            "--bowtie2db", str(db_dir), "--index", index, "--nproc", str(threads),
            "--unclassified_estimation", "-o", str(profile_path),
        ]
    else:
        if bowtie2out.is_file():
            bowtie2out.unlink()  # re-aligning; MetaPhlAn refuses to overwrite it
        command += ["--bowtie2out", str(bowtie2out)]
        if sam_out is not None:
            sam_out.parent.mkdir(parents=True, exist_ok=True)
            command += ["-s", str(sam_out)]

    reporter.info(
        f"  MetaPhlAn {tools.version} against {index} (extended catalogue"
        + (", with the marker SAM for strain typing)" if need_sam else ")")
    )
    started = time.monotonic()
    result = subprocess.run(
        command, capture_output=True, text=True, check=False, env=_env_for(tools.metaphlan)
    )
    elapsed = time.monotonic() - started
    if result.returncode != 0 or not profile_path.is_file():
        raise OpenBiotaError(
            f"MetaPhlAn 4 failed (exit {result.returncode}):\n{result.stderr[-3000:]}"
        )
    if need_sam and not (sam_out.is_file() and sam_out.stat().st_size > 0):  # type: ignore[union-attr]
        raise OpenBiotaError(
            f"MetaPhlAn 4 finished but wrote no marker SAM at {sam_out}; strain typing "
            "cannot run without it"
        )
    clades, unclassified, n_reads = parse_profile4(profile_path)
    stamp.write_text(
        json.dumps(
            {
                "cache_version": CACHE_VERSION,
                "inputs": fingerprint,
                "version": tools.version,
                "index": index,
                "elapsed_s": round(elapsed, 1),
                "command": [Path(command[0]).name, *command[1:]],
            }
        )
    )
    catalogue = ExtendedCatalogue(
        index=index, profiler_version=tools.version, clades=clades,
        unclassified_percent=unclassified, n_reads_processed=n_reads, elapsed_s=elapsed,
        cached=False, command=tuple(command),
    )
    reporter.ok(
        f"extended catalogue: {catalogue.n_sgbs} SGBs ({catalogue.n_unnamed} unnamed), "
        f"{catalogue.n_species} species, {unclassified:.1f}% unclassified, "
        f"{human_duration(elapsed)}"
    )
    return catalogue


__all__ = [
    "PINNED_INDEX4",
    "RENAMED",
    "ExtendedCatalogue",
    "Metaphlan4Tools",
    "SGBRow",
    "database4_ready",
    "locate_tools4",
    "parse_profile4",
    "run_metaphlan4",
]


# --------------------------------------------------------------------------- #
# Strain typing: marker consensus (sample2markers)
# --------------------------------------------------------------------------- #

#: Where `sample2markers` leaves its per-sample consensus, under the run's
#: ``strain/`` directory. `cli._strain_block` reads exactly this.
CONSENSUS_DIRNAME: Final = "consensus_markers"


def locate_sample2markers(tools: Metaphlan4Tools) -> str:
    """The `sample2markers.py` that ships beside the located MetaPhlAn 4."""
    candidate = Path(tools.metaphlan).parent / "sample2markers.py"
    if candidate.is_file():
        return str(candidate)
    found = shutil.which("sample2markers.py")
    if found:
        return found
    raise DependencyError(
        f"sample2markers.py not found beside {tools.metaphlan} or on PATH; it is part of "
        "the metaphlan>=4.1 package that strain typing needs"
    )


def run_sample2markers(
    tools: Metaphlan4Tools,
    *,
    sam: Path,
    strain_dir: Path,
    db_dir: Path,
    reporter: Reporter,
    threads: int = 8,
    index: str = PINNED_INDEX4,
) -> Path:
    """Build the dominant-population marker consensus for one sample.

    This is the second pass strain typing depends on. It used to be run by
    hand from a shell script after the pipeline, which is how one sample's
    report came to say "strain typing did not run" while the run itself
    reported success. It is a stage now, and a missing output is an error.
    """
    if shutil.which("samtools") is None:
        raise DependencyError(
            "samtools is not on PATH; sample2markers shells out to it "
            "(brew install samtools, or `make strain-system-deps`)"
        )
    if not (sam.is_file() and sam.stat().st_size > 0):
        raise OpenBiotaError(f"no marker SAM at {sam}; MetaPhlAn 4 must run with -s first")
    pkl = db_dir / f"{index}.pkl"
    if not pkl.is_file():
        raise DependencyError(f"MetaPhlAn 4 database pickle missing: {pkl}")
    out_dir = strain_dir / CONSENSUS_DIRNAME
    out_dir.mkdir(parents=True, exist_ok=True)
    stamp = strain_dir / "consensus_markers.json"
    existing = sorted(out_dir.glob("*.json.bz2"))
    if existing and stamp.is_file():
        meta = json.loads(stamp.read_text())
        if meta.get("cache_version") == CACHE_VERSION and meta.get("sam_size") == sam.stat().st_size:
            reporter.record(f"    strain markers: cached ({existing[0].name})")
            return existing[0]

    script = locate_sample2markers(tools)
    command = [
        script, "-i", str(sam), "-o", str(out_dir), "-d", str(pkl), "-n", str(threads),
    ]
    reporter.info(f"  sample2markers against {index}")
    started = time.monotonic()
    log = strain_dir / f"{sam.name.split('.')[0]}.s2m.log"
    result = subprocess.run(
        command, capture_output=True, text=True, check=False, env=_env_for(tools.metaphlan),
    )
    log.write_text((result.stdout or "") + (result.stderr or ""))
    elapsed = time.monotonic() - started
    produced = sorted(out_dir.glob("*.json.bz2"))
    if result.returncode != 0 or not produced:
        raise OpenBiotaError(
            f"sample2markers failed (exit {result.returncode}, log {log.name}):\n"
            f"{(result.stderr or result.stdout)[-2000:]}"
        )
    stamp.write_text(json.dumps({
        "cache_version": CACHE_VERSION, "sam_size": sam.stat().st_size,
        "elapsed_s": round(elapsed, 1), "command": [Path(command[0]).name, *command[1:]],
    }))
    reporter.ok(f"strain markers: {produced[0].name}, {human_duration(elapsed)}")
    return produced[0]
