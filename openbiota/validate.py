"""Validation: measure what this pipeline gets right and wrong, against truth.

The problem with validating a metagenomic screen is knowing the right answer.
This module solves that by constructing the sample itself.

1.  **Truth** comes from complete reference genomes. For each genome, its own
    annotated proteome is searched at *full length* against the panel
    reference sets, requiring a genuine ortholog — high identity across most
    of both proteins — plus the active-site residue criterion where a panel
    defines one. That is the published method's decision procedure, applied
    where it works properly: on whole proteins.

2.  **The test sample** is simulated from those same genomes at known cell
    abundances, so the true number of gene copies per 100 genomes is known
    exactly rather than estimated.

3.  **The measurement** is this pipeline, run on the simulated reads with no
    knowledge of the above.

Comparing 3 against 1 and 2 gives detection sensitivity, the false-positive
rate, and a calibration curve for the normalisation arithmetic. Four
experiments are run:

* ``balanced``  — 20 genomes at even abundance: sensitivity and specificity.
* ``negative``  — only genomes that truly lack the target genes: the
  false-positive rate, measured directly.
* ``calibration`` — known versus reported copies per 100 genomes.
* ``detection_limit`` — one carrier spiked down through 5%, 1% and 0.2% to
  find where detection fails.
"""

from __future__ import annotations

import argparse
import io
import json
import shutil
import statistics
import subprocess
import time
import zipfile
from collections.abc import Sequence
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Final

from openbiota.confidence import MIN_FRAGMENTS, MIN_MEAN_IDENTITY
from openbiota.errors import OpenBiotaError
from openbiota.logging_util import Reporter, human_bytes, human_duration
from openbiota.net import get
from openbiota.panels import PanelSet, load_panel_set
from openbiota.references import ReferenceDatabase, build_reference_database
from openbiota.search import SearchConfig, discover_sample, output_fields
from openbiota.seqio import iter_fasta, write_fasta
from openbiota.simulate import GenomeSpec, write_community

NCBI_DATASETS: Final = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"

#: Full-length ortholog criterion for ground truth. A protein counts as the
#: gene when it matches a curated reference at this identity across most of
#: both sequences. Well above the read-level floor, because at full length
#: there is no excuse for ambiguity.
TRUTH_MIN_IDENTITY: Final = 70.0
TRUTH_MIN_SUBJECT_COVERAGE: Final = 70.0
TRUTH_MIN_QUERY_COVERAGE: Final = 60.0


# --------------------------------------------------------------------------- #
# the validation panel of genomes
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class ValidationGenome:
    """One reference genome, with why it is in the set."""

    label: str
    species: str
    role: str  # "carrier" | "background"
    note: str = ""


#: Chosen so that every panel has both organisms that carry its gene and
#: organisms that do not. Species are named rather than pinned to accessions
#: so the set stays valid as assemblies are updated; the resolved accession is
#: recorded in the results for reproducibility.
VALIDATION_GENOMES: Final = (
    # --- known carriers, one or more panels each ---
    ValidationGenome("cscindens", "Clostridium scindens", "carrier", "bai operon reference organism"),
    ValidationGenome("chylemonae", "Clostridium hylemonae", "carrier", "bai operon"),
    ValidationGenome("phiranonis", "Peptacetobacter hiranonis", "carrier", "bai operon"),
    ValidationGenome("csporogenes", "Clostridium sporogenes", "carrier", "fld pathway, cutC"),
    ValidationGenome("cdifficile", "Clostridioides difficile", "carrier", "hpdBCA, p-cresol"),
    ValidationGenome("kpneumoniae", "Klebsiella pneumoniae", "carrier", "cutC, cntA"),
    ValidationGenome("ecoli", "Escherichia coli", "carrier", "tnaA, cntA/yeaW"),
    ValidationGenome("pmirabilis", "Proteus mirabilis", "carrier", "cutC"),
    ValidationGenome("ddesulfuricans", "Desulfovibrio desulfuricans", "carrier", "cutC"),
    ValidationGenome("fprausnitzii", "Faecalibacterium prausnitzii", "carrier", "but, butyrate"),
    ValidationGenome("rintestinalis", "Roseburia intestinalis", "carrier", "but, butyrate"),
    ValidationGenome("cbutyricum", "Clostridium butyricum", "carrier", "buk, butyrate"),
    ValidationGenome("elenta", "Eggerthella lenta", "carrier", "urdA, imidazole propionate"),
    ValidationGenome("spasteurianus", "Streptococcus pasteurianus", "carrier", "urdA — the source study's key organism"),
    ValidationGenome("lplantarum", "Lactiplantibacillus plantarum", "carrier", "urdA"),
    ValidationGenome("melsdenii", "Megasphaera elsdenii", "carrier", "propionate, acrylate route"),
    ValidationGenome("senterica", "Salmonella enterica", "carrier", "pduP, propanediol route"),
    ValidationGenome("bfragilis", "Bacteroides fragilis", "carrier", "bsh"),
    # --- background: abundant gut organisms not known for the target genes ---
    ValidationGenome("btheta", "Bacteroides thetaiotaomicron", "background", "abundant Bacteroidota"),
    ValidationGenome("pcopri", "Segatella copri", "background", "abundant Bacteroidota"),
    ValidationGenome("amuciniphila", "Akkermansia muciniphila", "background", "Verrucomicrobiota"),
    ValidationGenome("badolescentis", "Bifidobacterium adolescentis", "background", "Actinomycetota"),
    ValidationGenome("blongum", "Bifidobacterium longum", "background", "Actinomycetota"),
    ValidationGenome("sthermophilus", "Streptococcus thermophilus", "background", "Bacillota"),
)

#: Genomes used for the negative-control community. Deliberately excludes
#: every known carrier, so anything a panel reports here is a false positive
#: — subject to the truth check below, which has the final say.
NEGATIVE_LABELS: Final = (
    "btheta", "pcopri", "amuciniphila", "badolescentis", "blongum", "sthermophilus",
)

#: Carrier spiked into a background community for the detection-limit test.
#: urdA is the interesting case: it is the primary target and the one whose
#: specificity is hardest.
SPIKE_LABEL: Final = "spasteurianus"
SPIKE_FRACTIONS: Final = (0.05, 0.01, 0.002)


# --------------------------------------------------------------------------- #
# genome acquisition
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class ResolvedGenome:
    spec: ValidationGenome
    accession: str
    organism: str
    genome_fasta: Path
    protein_fasta: Path
    genome_size: int = 0
    n_proteins: int = 0


def resolve_accession(species: str) -> tuple[str, str]:
    """Pick the best complete RefSeq assembly for a species name."""
    quoted = species.replace(" ", "%20")
    for filters in (
        "filters.assembly_level=complete_genome&filters.assembly_source=refseq",
        "filters.assembly_source=refseq",
        "",
    ):
        url = f"{NCBI_DATASETS}/genome/taxon/{quoted}/dataset_report?page_size=20&{filters}"
        try:
            payload = json.loads(get(url, timeout=120).text())
        except (OpenBiotaError, json.JSONDecodeError):
            continue
        reports = payload.get("reports") or []
        refseq = [r for r in reports if str(r.get("accession", "")).startswith("GCF_")]
        for report in refseq or reports:
            accession = report.get("accession")
            if accession:
                return accession, report.get("organism", {}).get("organism_name", species)
    raise OpenBiotaError(
        f"no assembly found for {species!r} in NCBI Datasets. The species may have been "
        "renamed; update VALIDATION_GENOMES."
    )


def download_genome(
    accession: str, dest_dir: Path, *, refresh: bool = False
) -> tuple[Path, Path]:
    """Fetch genomic and protein FASTA for one assembly."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    genome_path = dest_dir / f"{accession}_genomic.fna"
    protein_path = dest_dir / f"{accession}_protein.faa"
    if genome_path.is_file() and protein_path.is_file() and not refresh:
        return genome_path, protein_path

    url = (
        f"{NCBI_DATASETS}/genome/accession/{accession}/download"
        "?include_annotation_type=GENOME_FASTA&include_annotation_type=PROT_FASTA"
    )
    body = get(url, timeout=900).body
    try:
        archive = zipfile.ZipFile(io.BytesIO(body))
    except zipfile.BadZipFile as exc:
        raise OpenBiotaError(f"{accession}: NCBI returned a non-zip response ({len(body)} bytes)") from exc

    genomic = [n for n in archive.namelist() if n.endswith(("_genomic.fna", ".fna"))]
    proteins = [n for n in archive.namelist() if n.endswith(("protein.faa", ".faa"))]
    if not genomic:
        raise OpenBiotaError(f"{accession}: no genomic FASTA in the NCBI archive")
    if not proteins:
        raise OpenBiotaError(
            f"{accession}: no protein FASTA in the NCBI archive — the assembly is unannotated "
            "and cannot provide ground truth"
        )

    genome_path.write_bytes(archive.read(genomic[0]))
    protein_path.write_bytes(archive.read(proteins[0]))
    return genome_path, protein_path


def acquire_genomes(
    genomes: Sequence[ValidationGenome],
    dest_dir: Path,
    reporter: Reporter,
    *,
    refresh: bool = False,
) -> list[ResolvedGenome]:
    cache = dest_dir / "accessions.json"
    resolved_cache: dict[str, list[str]] = {}
    if cache.is_file() and not refresh:
        resolved_cache = json.loads(cache.read_text(encoding="utf-8"))

    out: list[ResolvedGenome] = []
    for spec in genomes:
        if spec.label in resolved_cache:
            accession, organism = resolved_cache[spec.label]
        else:
            accession, organism = resolve_accession(spec.species)
            resolved_cache[spec.label] = [accession, organism]
            cache.parent.mkdir(parents=True, exist_ok=True)
            cache.write_text(json.dumps(resolved_cache, indent=2), encoding="utf-8")
            reporter.record(f"    {spec.label}: {accession} ({organism})")

        genome_fasta, protein_fasta = download_genome(
            accession, dest_dir / spec.label, refresh=refresh
        )
        size = sum(len(s) for _, s in iter_fasta(genome_fasta))
        n_proteins = sum(1 for _ in iter_fasta(protein_fasta))
        out.append(
            ResolvedGenome(
                spec=spec,
                accession=accession,
                organism=organism,
                genome_fasta=genome_fasta,
                protein_fasta=protein_fasta,
                genome_size=size,
                n_proteins=n_proteins,
            )
        )
        reporter.info(
            f"  {spec.label:<16} {accession:<18} {size / 1e6:>5.2f} Mb, "
            f"{n_proteins:>5,} proteins  {organism[:40]}"
        )
    return out


# --------------------------------------------------------------------------- #
# ground truth from full-length proteomes
# --------------------------------------------------------------------------- #


@dataclass(frozen=True, slots=True)
class TruthCall:
    """Gene presence in one genome, at full length.

    Two counts, because two different things are being validated:

    * ``copies`` applies the active-site residue criterion where a panel
      defines one. This is the published method's answer, and it is what the
      pipeline's residue-consistent subset should be compared against.
    * ``copies_any_homolog`` counts full-length orthologs without the residue
      criterion. This is what the pipeline's *headline* figure should be
      compared against, since the headline does not apply that filter either.

    Comparing the headline against the residue-filtered truth would be
    comparing two different quantities and would misattribute the difference
    to inaccuracy.
    """

    genome: str
    entry_key: str
    present: bool
    copies: int
    present_any_homolog: bool
    copies_any_homolog: int
    best_identity: float
    best_subject_coverage: float
    residue_ok: bool | None
    accession: str = ""

    def to_json(self) -> dict[str, Any]:
        return asdict(self)


def determine_truth(
    genomes: Sequence[ResolvedGenome],
    database: ReferenceDatabase,
    panel_set: PanelSet,
    work_dir: Path,
    reporter: Reporter,
    *,
    diamond: str = "diamond",
    threads: int = 1,
) -> dict[tuple[str, str], TruthCall]:
    """Call gene presence per genome from its proteome, at full length.

    This is the reference method: whole predicted proteins, a strict ortholog
    threshold, and the active-site residue check applied where the panel
    defines one — exactly where that check is decisive.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    target_keys = {e.key for p in panel_set.panels for e in p.targets}

    # Build a DIAMOND database of just the target references.
    target_fasta = work_dir / "targets.faa"
    records: list[tuple[str, str]] = []
    seq_by_idx = {int(h.split("~")[0]): s for h, s in iter_fasta(database.fasta_path)}
    for reference in database.sequences:
        if reference.entry_key in target_keys:
            records.append((reference.sseqid, seq_by_idx[reference.idx]))
    write_fasta(target_fasta, records)

    target_db = work_dir / "targets"
    _run([_diamond(diamond), "makedb", "--in", str(target_fasta), "--db", str(target_db),
          "--threads", str(threads), "--quiet"], reporter)

    truth: dict[tuple[str, str], TruthCall] = {}
    anchor_positions = {
        idx: pos
        for mapping in database.anchors.values()
        for idx, pos in mapping.positions.items()
    }
    residue_panels = {
        p.name: p.residue_check for p in panel_set.panels if p.residue_check and p.residue_check.enabled
    }

    for genome in genomes:
        out_tsv = work_dir / f"{genome.spec.label}_truth.tsv"
        _run(
            [
                _diamond(diamond), "blastp",
                "--db", str(target_db),
                "--query", str(genome.protein_fasta),
                "--out", str(out_tsv),
                "--outfmt", "6", "qseqid", "sseqid", "pident", "qlen", "slen",
                "length", "sstart", "send", "qseq_gapped", "sseq_gapped",
                "--max-target-seqs", "1",
                "--index-chunks", "1",
                "--very-sensitive",
                "--evalue", "1e-10",
                "--threads", str(threads),
                "--quiet",
            ],
            reporter,
        )

        # entry_key -> list of (identity, subject_coverage, residue_ok)
        hits: dict[str, list[tuple[float, float, bool | None]]] = {}
        with out_tsv.open("r", encoding="utf-8") as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 8:
                    continue
                try:
                    pident = float(parts[2])
                    qlen, slen, aln = int(parts[3]), int(parts[4]), int(parts[5])
                    sstart, send = int(parts[6]), int(parts[7])
                except ValueError:
                    continue
                reference = database.by_sseqid(parts[1])
                if reference is None or reference.entry_key not in target_keys:
                    continue
                s_cov = aln / slen * 100.0 if slen else 0.0
                q_cov = aln / qlen * 100.0 if qlen else 0.0
                if (
                    pident < TRUTH_MIN_IDENTITY
                    or s_cov < TRUTH_MIN_SUBJECT_COVERAGE
                    or q_cov < TRUTH_MIN_QUERY_COVERAGE
                ):
                    continue

                residue_ok: bool | None = None
                check = residue_panels.get(reference.panel)
                if check is not None and reference.entry_id in check.applies_to:
                    position = anchor_positions.get(reference.idx)
                    residue_ok = False
                    if position is not None and len(parts) >= 10:
                        from openbiota.residues import residue_at_subject_position

                        residue = residue_at_subject_position(
                            parts[8], parts[9], min(sstart, send), position
                        )
                        residue_ok = bool(
                            residue and residue.upper() in check.accepted_residues
                        )
                hits.setdefault(reference.entry_key, []).append((pident, s_cov, residue_ok))

        for key in target_keys:
            found = hits.get(key, [])
            # A panel with a residue check requires it to pass at full length —
            # this is the criterion the published method applies.
            qualifying = [h for h in found if h[2] is not False]
            best = max(found, key=lambda h: h[0]) if found else (0.0, 0.0, None)
            truth[(genome.spec.label, key)] = TruthCall(
                genome=genome.spec.label,
                entry_key=key,
                present=bool(qualifying),
                copies=len(qualifying),
                present_any_homolog=bool(found),
                copies_any_homolog=len(found),
                best_identity=best[0],
                best_subject_coverage=best[1],
                residue_ok=best[2],
                accession=genome.accession,
            )

    n_present = sum(1 for t in truth.values() if t.present)
    n_any = sum(1 for t in truth.values() if t.present_any_homolog)
    reporter.info(
        f"  ground truth: {n_present} gene/genome pairs present out of {len(truth)} tested "
        f"({n_any} counting full-length homologs that fail an active-site criterion)"
    )
    return truth


def _diamond(diamond: str) -> str:
    exe = shutil.which(diamond)
    if exe is None:
        raise OpenBiotaError(f"{diamond!r} not found on PATH")
    return exe


def _run(cmd: list[str], reporter: Reporter) -> None:
    reporter.record(f"    $ {' '.join(cmd)}")
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, timeout=7200)  # noqa: S603
    except subprocess.CalledProcessError as exc:
        raise OpenBiotaError(
            f"command failed: {' '.join(cmd)}\n  stderr: {(exc.stderr or '').strip()[-1500:]}"
        ) from exc


# --------------------------------------------------------------------------- #
# expected values
# --------------------------------------------------------------------------- #


def expected_copies_per_100(
    community: Sequence[GenomeSpec],
    truth: dict[tuple[str, str], TruthCall],
    entry_key: str,
    *,
    residue_filtered: bool = True,
) -> float:
    """True copies per 100 bacterial genomes for one gene in one community.

    rpoB normalisation counts *genome equivalents*, so the weighting is by cell
    abundance, not by DNA mass. ``residue_filtered`` selects which of the two
    truth counts to use — see :class:`TruthCall`.
    """
    total_cells = sum(g.abundance for g in community)
    if total_cells <= 0:
        return 0.0
    copies = 0.0
    for member in community:
        call = truth.get((member.name, entry_key))
        if call is not None:
            copies += member.abundance * (
                call.copies if residue_filtered else call.copies_any_homolog
            )
    return copies / total_cells * 100.0


# --------------------------------------------------------------------------- #
# experiments
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class GeneOutcome:
    entry_key: str
    panel: str
    gene: str
    truly_present: bool
    true_copies_per_100: float
    detected: bool
    reported_fragments: int
    reported_copies_per_100: float | None
    mean_identity: float | None
    classification: str = ""
    true_copies_residue_filtered: float = 0.0
    reported_copies_residue_consistent: float | None = None

    def to_json(self) -> dict[str, Any]:
        return asdict(self)


@dataclass(slots=True)
class ExperimentResult:
    name: str
    description: str
    n_read_pairs: int
    community: list[dict[str, Any]]
    rpob_fragments: int
    outcomes: list[GeneOutcome] = field(default_factory=list)

    def confusion(self) -> dict[str, int]:
        counts = {"TP": 0, "FP": 0, "FN": 0, "TN": 0}
        for outcome in self.outcomes:
            counts[outcome.classification] = counts.get(outcome.classification, 0) + 1
        return counts

    def to_json(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "description": self.description,
            "n_read_pairs": self.n_read_pairs,
            "community": self.community,
            "rpob_fragments": self.rpob_fragments,
            "confusion": self.confusion(),
            "genes": [o.to_json() for o in self.outcomes],
        }


def classify(truly_present: bool, detected: bool) -> str:
    if truly_present and detected:
        return "TP"
    if truly_present and not detected:
        return "FN"
    if not truly_present and detected:
        return "FP"
    return "TN"


def _metrics(counts: dict[str, int]) -> dict[str, float | None]:
    tp, fp, fn, tn = counts.get("TP", 0), counts.get("FP", 0), counts.get("FN", 0), counts.get("TN", 0)
    def ratio(num: int, den: int) -> float | None:
        return num / den if den else None
    return {
        "sensitivity": ratio(tp, tp + fn),
        "specificity": ratio(tn, tn + fp),
        "precision": ratio(tp, tp + fp),
        "false_positive_rate": ratio(fp, fp + tn),
        "false_discovery_rate": ratio(fp, tp + fp),
        "accuracy": ratio(tp + tn, tp + tn + fp + fn),
    }


def _linear_fit(pairs: Sequence[tuple[float, float]]) -> dict[str, float | None]:
    """Least-squares fit of reported against true, plus R²."""
    usable = [(x, y) for x, y in pairs if x > 0 or y > 0]
    if len(usable) < 3:
        return {"slope": None, "intercept": None, "r_squared": None, "n": len(usable)}
    xs = [x for x, _ in usable]
    ys = [y for _, y in usable]
    mean_x, mean_y = statistics.fmean(xs), statistics.fmean(ys)
    sxx = sum((x - mean_x) ** 2 for x in xs)
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in usable)
    if sxx == 0:
        return {"slope": None, "intercept": None, "r_squared": None, "n": len(usable)}
    slope = sxy / sxx
    intercept = mean_y - slope * mean_x
    ss_tot = sum((y - mean_y) ** 2 for y in ys)
    ss_res = sum((y - (slope * x + intercept)) ** 2 for x, y in usable)
    r2 = 1 - ss_res / ss_tot if ss_tot else None
    return {
        "slope": slope,
        "intercept": intercept,
        "r_squared": r2,
        "n": len(usable),
    }


# --------------------------------------------------------------------------- #
# orchestration
# --------------------------------------------------------------------------- #


def _screen_community(
    *,
    community: Sequence[GenomeSpec],
    name: str,
    n_pairs: int,
    seed: int,
    error_rate: float,
    work_dir: Path,
    database: ReferenceDatabase,
    panel_set: PanelSet,
    config: SearchConfig,
    fields: Sequence[str],
    reporter: Reporter,
) -> Any:
    """Simulate a community, run the pipeline over it, return the tally."""
    from openbiota.cli import screen_one

    fastq_dir = work_dir / "simulated"
    started = time.monotonic()
    r1, r2, written = write_community(
        community,
        out_prefix=fastq_dir / name,
        n_pairs=n_pairs,
        seed=seed,
        error_rate=error_rate,
    )
    reporter.info(
        f"  {name}: simulated {written:,} read pairs from {len(community)} genomes "
        f"({human_bytes(r1.stat().st_size + r2.stat().st_size)}) in "
        f"{human_duration(time.monotonic() - started)}"
    )

    sample_input = discover_sample(r1=r1, r2=r2, sample=name)
    return screen_one(
        database=database,
        panel_set=panel_set,
        sample_input=sample_input,
        out_dir=work_dir / "results" / name,
        config=config,
        fields=fields,
        reporter=reporter,
    )


def _outcomes_for(
    result: Any,
    community: Sequence[GenomeSpec],
    truth: dict[tuple[str, str], TruthCall],
) -> list[GeneOutcome]:
    outcomes: list[GeneOutcome] = []
    for panel in result.panels:
        has_residue_check = any(t.residue is not None for t in panel.targets)
        for gene in panel.targets:
            # The headline figure does not apply an active-site filter, so it
            # is compared against homolog-level truth. The residue-consistent
            # subset is compared against residue-filtered truth.
            true_headline = expected_copies_per_100(
                community, truth, gene.key, residue_filtered=not has_residue_check
            )
            true_strict = expected_copies_per_100(
                community, truth, gene.key, residue_filtered=True
            )
            truly_present = true_headline > 0
            detected = gene.fragments > 0
            outcomes.append(
                GeneOutcome(
                    entry_key=gene.key,
                    panel=panel.panel.name,
                    gene=gene.gene or gene.entry_id,
                    truly_present=truly_present,
                    true_copies_per_100=round(true_headline, 4),
                    detected=detected,
                    reported_fragments=gene.fragments,
                    reported_copies_per_100=(
                        None
                        if gene.copies_per_100_genomes is None
                        else round(gene.copies_per_100_genomes, 4)
                    ),
                    mean_identity=(
                        None if gene.mean_identity is None else round(gene.mean_identity, 1)
                    ),
                    classification=classify(truly_present, detected),
                    true_copies_residue_filtered=round(true_strict, 4),
                    reported_copies_residue_consistent=(
                        None
                        if gene.copies_per_100_genomes_residue_consistent is None
                        else round(gene.copies_per_100_genomes_residue_consistent, 4)
                    ),
                )
            )
    return outcomes


def run_validation(args: argparse.Namespace) -> int:
    """Execute the four validation experiments and write the results."""
    reporter = Reporter(verbose=not args.quiet)
    work_dir = Path(args.work_dir)
    panel_set = load_panel_set(args.panels_dir)

    with reporter.stage("reference database"):
        database = build_reference_database(
            panel_set,
            [panel_set.normalizer, *(e for p in panel_set.panels for e in p.entries)],
            refs_dir=args.refs_dir,
            reporter=reporter,
            diamond=args.diamond,
            threads=args.threads,
        )

    with reporter.stage(f"acquiring {len(VALIDATION_GENOMES)} reference genomes"):
        genomes = acquire_genomes(
            VALIDATION_GENOMES, work_dir / "genomes", reporter, refresh=args.refresh
        )
    by_label = {g.spec.label: g for g in genomes}

    with reporter.stage("ground truth from full-length proteomes"):
        truth = determine_truth(
            genomes,
            database,
            panel_set,
            work_dir / "truth",
            reporter,
            diamond=args.diamond,
            threads=args.threads,
        )

    fields = output_fields(with_residues=bool(database.anchors))
    config = SearchConfig(
        diamond=args.diamond, threads=args.threads, block_size=4.0, index_chunks=1
    )
    experiments: list[ExperimentResult] = []

    def community_json(community: Sequence[GenomeSpec]) -> list[dict[str, Any]]:
        total = sum(g.abundance for g in community)
        return [
            {
                "label": g.name,
                "organism": by_label[g.name].organism,
                "accession": by_label[g.name].accession,
                "cell_fraction": round(g.abundance / total, 6),
                "genome_size": by_label[g.name].genome_size,
            }
            for g in community
        ]

    # ---- experiment 1: balanced community ------------------------------------ #
    balanced = [
        GenomeSpec(name=g.spec.label, path=g.genome_fasta, abundance=1.0) for g in genomes
    ]
    with reporter.stage("experiment 1/4 — balanced community"):
        result = _screen_community(
            community=balanced,
            name="balanced",
            n_pairs=args.reads,
            seed=args.seed,
            error_rate=args.error_rate,
            work_dir=work_dir,
            database=database,
            panel_set=panel_set,
            config=config,
            fields=fields,
            reporter=reporter,
        )
        experiments.append(
            ExperimentResult(
                name="balanced",
                description=(
                    f"All {len(balanced)} validation genomes at equal cell abundance. "
                    "Measures detection sensitivity and specificity when every gene that "
                    "is present is present at a substantial fraction of the community."
                ),
                n_read_pairs=args.reads,
                community=community_json(balanced),
                rpob_fragments=result.normalizer.fragments,
                outcomes=_outcomes_for(result, balanced, truth),
            )
        )

    # ---- experiment 2: negative control -------------------------------------- #
    negative = [
        GenomeSpec(name=by_label[label].spec.label, path=by_label[label].genome_fasta, abundance=1.0)
        for label in NEGATIVE_LABELS
        if label in by_label
    ]
    with reporter.stage("experiment 2/4 — negative control"):
        result = _screen_community(
            community=negative,
            name="negative",
            n_pairs=args.reads,
            seed=args.seed + 1,
            error_rate=args.error_rate,
            work_dir=work_dir,
            database=database,
            panel_set=panel_set,
            config=config,
            fields=fields,
            reporter=reporter,
        )
        experiments.append(
            ExperimentResult(
                name="negative",
                description=(
                    f"{len(negative)} abundant gut organisms with no known carriage of the "
                    "target genes. Gene presence is verified from each proteome, so anything "
                    "reported for an absent gene is a false positive. This is the experiment "
                    "that measures the false-positive rate directly."
                ),
                n_read_pairs=args.reads,
                community=community_json(negative),
                rpob_fragments=result.normalizer.fragments,
                outcomes=_outcomes_for(result, negative, truth),
            )
        )

    # ---- experiment 3: detection limit --------------------------------------- #
    spike_results: list[dict[str, Any]] = []
    if SPIKE_LABEL in by_label:
        for index, fraction in enumerate(SPIKE_FRACTIONS):
            background = [
                GenomeSpec(
                    name=by_label[label].spec.label,
                    path=by_label[label].genome_fasta,
                    abundance=(1.0 - fraction) / len(NEGATIVE_LABELS),
                )
                for label in NEGATIVE_LABELS
                if label in by_label
            ]
            community = [
                *background,
                GenomeSpec(
                    name=SPIKE_LABEL,
                    path=by_label[SPIKE_LABEL].genome_fasta,
                    abundance=fraction,
                ),
            ]
            name = f"spike_{fraction:g}".replace(".", "p")
            with reporter.stage(
                f"experiment 3/4 — detection limit, carrier at {fraction:.1%} "
                f"({index + 1}/{len(SPIKE_FRACTIONS)})"
            ):
                result = _screen_community(
                    community=community,
                    name=name,
                    n_pairs=args.reads,
                    seed=args.seed + 10 + index,
                    error_rate=args.error_rate,
                    work_dir=work_dir,
                    database=database,
                    panel_set=panel_set,
                    config=config,
                    fields=fields,
                    reporter=reporter,
                )
            outcomes = _outcomes_for(result, community, truth)
            experiments.append(
                ExperimentResult(
                    name=name,
                    description=(
                        f"{by_label[SPIKE_LABEL].organism} spiked at {fraction:.2%} of cells "
                        "into a background of non-carriers. Finds the abundance at which "
                        "detection fails."
                    ),
                    n_read_pairs=args.reads,
                    community=community_json(community),
                    rpob_fragments=result.normalizer.fragments,
                    outcomes=outcomes,
                )
            )
            urda = next((o for o in outcomes if o.entry_key == "urda:URDA"), None)
            if urda is not None:
                spike_results.append(
                    {
                        "carrier_cell_fraction": fraction,
                        "true_copies_per_100": urda.true_copies_per_100,
                        "detected": urda.detected,
                        "reported_fragments": urda.reported_fragments,
                        "reported_copies_per_100": urda.reported_copies_per_100,
                    }
                )

    # ---- experiment 4: calibration ------------------------------------------- #
    calibration_points = [
        (o.true_copies_per_100, o.reported_copies_per_100)
        for experiment in experiments
        for o in experiment.outcomes
        if o.reported_copies_per_100 is not None
    ]
    calibration = _linear_fit(calibration_points)

    # Per-gene accuracy. The systematic component matters more than the scatter:
    # a bias that is constant per gene cancels when the same gene is compared
    # across samples, which is why results are reported as percentiles.
    per_gene_ratio: dict[str, list[float]] = {}
    for experiment in experiments:
        for o in experiment.outcomes:
            if (
                o.reported_copies_per_100 is not None
                and o.true_copies_per_100 > 0
                and o.reported_fragments >= MIN_FRAGMENTS
                and (o.mean_identity or 0) >= MIN_MEAN_IDENTITY
            ):
                per_gene_ratio.setdefault(o.entry_key, []).append(
                    o.reported_copies_per_100 / o.true_copies_per_100
                )
    per_gene_calibration = {
        key: {
            "n": len(values),
            "median_reported_over_true": round(statistics.median(values), 3),
            "min": round(min(values), 3),
            "max": round(max(values), 3),
        }
        for key, values in sorted(per_gene_ratio.items())
    }

    # Does the residue-consistent subset track residue-filtered truth better
    # than the headline tracks it? This is the test of whether that filter works.
    residue_pairs_headline: list[tuple[float, float]] = []
    residue_pairs_strict: list[tuple[float, float]] = []
    for experiment in experiments:
        for o in experiment.outcomes:
            if o.reported_copies_residue_consistent is None:
                continue
            if o.reported_copies_per_100 is not None:
                residue_pairs_headline.append(
                    (o.true_copies_residue_filtered, o.reported_copies_per_100)
                )
            residue_pairs_strict.append(
                (o.true_copies_residue_filtered, o.reported_copies_residue_consistent)
            )
    residue_effect = {
        "headline_vs_residue_filtered_truth": _linear_fit(residue_pairs_headline),
        "residue_consistent_vs_residue_filtered_truth": _linear_fit(residue_pairs_strict),
    }

    # ---- aggregate ----------------------------------------------------------- #
    # Two operating points. "permissive" counts any fragment as a detection;
    # "confirmed" applies the fragment-count and identity criteria that
    # openbiota.confidence uses in production. Reporting both is what makes the
    # confirmed thresholds defensible rather than arbitrary.
    def confirmed_call(outcome: GeneOutcome) -> bool:
        return (
            outcome.reported_fragments >= MIN_FRAGMENTS
            and (outcome.mean_identity or 0.0) >= MIN_MEAN_IDENTITY
        )

    overall: dict[str, int] = {"TP": 0, "FP": 0, "FN": 0, "TN": 0}
    overall_confirmed: dict[str, int] = {"TP": 0, "FP": 0, "FN": 0, "TN": 0}
    per_panel: dict[str, dict[str, int]] = {}
    per_panel_confirmed: dict[str, dict[str, int]] = {}
    for experiment in experiments:
        for outcome in experiment.outcomes:
            overall[outcome.classification] += 1
            per_panel.setdefault(outcome.panel, {"TP": 0, "FP": 0, "FN": 0, "TN": 0})[
                outcome.classification
            ] += 1

            tier = classify(outcome.truly_present, confirmed_call(outcome))
            overall_confirmed[tier] += 1
            per_panel_confirmed.setdefault(
                outcome.panel, {"TP": 0, "FP": 0, "FN": 0, "TN": 0}
            )[tier] += 1

    payload = {
        "tool": "openbiota",
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "reference_fingerprint": database.fingerprint,
        "diamond_version": database.diamond_version.splitlines()[0],
        "read_pairs_per_experiment": args.reads,
        "simulated_error_rate": args.error_rate,
        "seed": args.seed,
        "ground_truth_criterion": {
            "method": (
                "full-length blastp of each genome's annotated proteome against the panel "
                "reference sets, plus the active-site residue criterion where a panel "
                "defines one"
            ),
            "min_identity_percent": TRUTH_MIN_IDENTITY,
            "min_subject_coverage_percent": TRUTH_MIN_SUBJECT_COVERAGE,
            "min_query_coverage_percent": TRUTH_MIN_QUERY_COVERAGE,
        },
        "genomes": [
            {
                "label": g.spec.label,
                "species": g.spec.species,
                "organism": g.organism,
                "accession": g.accession,
                "role": g.spec.role,
                "note": g.spec.note,
                "genome_size": g.genome_size,
                "n_proteins": g.n_proteins,
            }
            for g in genomes
        ],
        "truth_calls": [t.to_json() for t in truth.values() if t.present],
        "experiments": [e.to_json() for e in experiments],
        "operating_points": {
            "permissive": {
                "criterion": "any matching fragment counts as a detection",
                "confusion": overall,
                "metrics": _metrics(overall),
            },
            "confirmed": {
                "criterion": (
                    f">= {MIN_FRAGMENTS} matching fragments and >= "
                    f"{MIN_MEAN_IDENTITY:.0f}% mean translated identity — the tier openbiota "
                    "reports as 'confirmed'"
                ),
                "confusion": overall_confirmed,
                "metrics": _metrics(overall_confirmed),
            },
        },
        "overall_confusion": overall_confirmed,
        "overall_metrics": _metrics(overall_confirmed),
        "per_panel_metrics": {
            k: _metrics(per_panel_confirmed.get(k, v)) | {"confusion": per_panel_confirmed.get(k, v)}
            for k, v in sorted(per_panel.items())
        },
        "per_panel_metrics_permissive": {
            k: _metrics(v) | {"confusion": v} for k, v in sorted(per_panel.items())
        },
        "calibration": calibration,
        "per_gene_calibration": per_gene_calibration,
        "residue_filter_effect": residue_effect,
        "detection_limit": spike_results,
    }

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    _print_validation_summary(payload)
    reporter.ok(f"validation results written to {out_path}")
    return 0


def _print_validation_summary(payload: dict[str, Any]) -> None:
    metrics = payload["overall_metrics"]
    counts = payload["overall_confusion"]

    def pct(value: float | None) -> str:
        return "n/a" if value is None else f"{value:.1%}"

    print()
    print("=" * 78)
    print(" VALIDATION RESULTS — measured against known ground truth")
    print("=" * 78)
    print()
    print(f"  genomes used            {len(payload['genomes'])}")
    print(f"  read pairs / experiment {payload['read_pairs_per_experiment']:,}")
    print(f"  experiments             {len(payload['experiments'])}")
    print()
    print("  OPERATING POINTS")
    print(f"    {'':<26}{'sens':>9}{'spec':>9}{'prec':>9}{'FP':>6}{'FN':>6}")
    print("    " + "-" * 65)
    for name, block in payload["operating_points"].items():
        m, c = block["metrics"], block["confusion"]
        print(
            f"    {name:<26}{pct(m['sensitivity']):>9}{pct(m['specificity']):>9}"
            f"{pct(m['precision']):>9}{c['FP']:>6}{c['FN']:>6}"
        )
    print()
    for name, block in payload["operating_points"].items():
        print(f"    {name}: {block['criterion']}")
    print()
    print("  AT THE CONFIRMED TIER (what openbiota reports)")
    print(f"    true positives        {counts['TP']}")
    print(f"    true negatives        {counts['TN']}")
    print(f"    false positives       {counts['FP']}")
    print(f"    false negatives       {counts['FN']}")
    print(f"    sensitivity           {pct(metrics['sensitivity'])}")
    print(f"    specificity           {pct(metrics['specificity'])}")
    print(f"    precision             {pct(metrics['precision'])}")
    print(f"    accuracy              {pct(metrics['accuracy'])}")
    print()
    calibration = payload["calibration"]
    if calibration.get("r_squared") is not None:
        print("  QUANTITATIVE CALIBRATION  (reported vs true copies per 100 genomes)")
        print(f"    slope                 {calibration['slope']:.3f}   (1.0 = perfect)")
        print(f"    intercept             {calibration['intercept']:.3f}")
        print(f"    R-squared             {calibration['r_squared']:.3f}")
        print(f"    points                {calibration['n']}")
        print()
    print(f"  {'panel':<14}{'sens':>8}{'spec':>8}{'prec':>8}{'FPR':>8}   confusion")
    print("  " + "-" * 66)
    for name, m in payload["per_panel_metrics"].items():
        c = m["confusion"]
        print(
            f"  {name:<14}{pct(m['sensitivity']):>8}{pct(m['specificity']):>8}"
            f"{pct(m['precision']):>8}{pct(m['false_positive_rate']):>8}"
            f"   TP={c['TP']} TN={c['TN']} FP={c['FP']} FN={c['FN']}"
        )
    per_gene = payload.get("per_gene_calibration") or {}
    if per_gene:
        print("  PER-GENE QUANTITATIVE BIAS  (reported / true, confirmed calls)")
        print(f"    {'gene':<20}{'median':>9}{'min':>8}{'max':>8}{'n':>5}")
        for key, stats in per_gene.items():
            print(
                f"    {key:<20}{stats['median_reported_over_true']:>9.2f}"
                f"{stats['min']:>8.2f}{stats['max']:>8.2f}{stats['n']:>5}"
            )
        print()

    if payload["detection_limit"]:
        print()
        print("  DETECTION LIMIT  (urdA carrier spiked into non-carrier background)")
        print(f"    {'cell fraction':>15}{'true c/100':>13}{'detected':>11}{'fragments':>11}")
        for row in payload["detection_limit"]:
            print(
                f"    {row['carrier_cell_fraction']:>14.2%}"
                f"{row['true_copies_per_100']:>13.2f}"
                f"{('YES' if row['detected'] else 'no'):>11}"
                f"{row['reported_fragments']:>11,}"
            )
    print()
