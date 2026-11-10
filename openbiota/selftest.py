"""Does every detection lane actually run on this machine?

    openbiota selftest                 # every detection lane, about ten minutes
    openbiota selftest --lanes singlem_globdb motus4

The installer used to call a lane installed when its files were on disk.
That is not the same thing as working: SingleM's files were all there on a
fresh Ubuntu machine and its first real run died on ``hmmsearch: command
not found``, fifteen minutes into somebody's report. A tool is installed
when it has run, so this runs them.

Reads are simulated from three genomes that travel with the repository
(the probiotic strain references in ``refs/panels/probiotics``), a few
times over each genome, and pushed through each lane by the pipeline's own
lane code - the same function ``openbiota run`` calls, with the same
environment. A lane passes when it runs to completion and finds what the
reads were drawn from; a lane that exits cleanly having seen nothing, or
something else, in reads from three well-catalogued gut species has not
worked either.

Each lane is run in a process of its own, which is what lets the line it
prints end with the most memory any one of its programs held: the figure
the documentation's memory requirement is made of.

"Finds" is judged by genus. The six catalogues agree that the reads are
*Bifidobacterium* and *Lacticaseibacillus* and on little else: the same
*Bifidobacterium animalis* genome is ``Bifidobacterium_animalis`` to
MetaPhlAn, ``Bifidobacterium MOTU40_005954`` to GlobDB and an unnamed mOTU
to mOTUs, and Kraken adds a hundred low-count species of its own.
"""

from __future__ import annotations

import gzip
import os
import shutil
import sys
import tempfile
import time
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Final

from openbiota.logging_util import Reporter, Style, human_duration, plain_reason

REPO: Final = Path(__file__).resolve().parents[1]

#: The genomes the reads are drawn from: (file under refs/panels/probiotics,
#: read-name tag, what it is). The tag goes into every read name, so it has no
#: spaces: a name is cut at the first space, and reads that then share a name
#: are counted once by the tools that count read pairs.
GENOMES: Final = (
    ("probiotic.bb12.fasta", "bb12", "Bifidobacterium animalis BB-12"),
    ("probiotic.blongum_35624.fasta", "blongum35624", "Bifidobacterium longum 35624"),
    ("probiotic.lgg_fm179322.fasta", "lgg", "Lacticaseibacillus rhamnosus GG"),
)

#: The genera those genomes belong to, in every catalogue's taxonomy.
EXPECTED_GENERA: Final = ("Bifidobacterium", "Lacticaseibacillus")

#: Read pairs simulated. About five-fold coverage of the three genomes:
#: enough for every lane's detection floor, small enough that the time goes
#: to loading each lane's database, which is the thing being tested.
N_PAIRS: Final = 120_000

#: lane id -> the name a person would recognise it by.
LANES: Final[dict[str, str]] = {
    "metaphlan_jan26": "MetaPhlAn 4 (Jan26 catalogue)",
    "sylph_globdb": "sylph (GlobDB)",
    "motus4": "mOTUs 4",
    "kraken_uhgg": "Kraken 2 + Bracken (UHGG)",
    "singlem_globdb": "SingleM (GlobDB)",
    "genome_profile": "sylph (GTDB whole-genome profile)",
}


@dataclass(frozen=True)
class LaneOutcome:
    lane: str
    label: str
    ok: bool
    seconds: float
    reason: str = ""
    #: The expected genera the lane named.
    found: tuple[str, ...] = ()
    #: The most memory any one of the lane's programs held, in GB; None when not measured.
    peak_gb: float | None = None


NOT_FOUND: Final = ("ran, but did not find the test organisms (" + ", ".join(EXPECTED_GENERA)
                    + ") in reads drawn from them")


def genera_found(lineages: Iterable[str]) -> tuple[str, ...]:
    """Which expected genera appear in these lineage strings, whatever the catalogue's spelling."""
    text = " ".join(lineages).replace("_", " ")
    return tuple(genus for genus in EXPECTED_GENERA if genus in text)


def simulate_reads(work: Path, *, n_pairs: int = N_PAIRS, seed: int = 20261110,
                   genomes_dir: Path | None = None) -> tuple[Path, Path]:
    """Paired, gzipped FASTQ drawn from the bundled genomes. Raises if they are absent.

    Gzipped because that is what a real run hands the lanes: the reads that
    come out of quality control and human-read removal.
    """
    from openbiota.simulate import GenomeSpec, write_community  # noqa: PLC0415

    folder = genomes_dir or (REPO / "refs" / "panels" / "probiotics")
    specs = [GenomeSpec(tag, folder / filename, 1.0) for filename, tag, _ in GENOMES]
    absent = [s.path.name for s in specs if not s.path.is_file()]
    if absent:
        raise FileNotFoundError(
            f"the self-test genomes are not in {folder} (" + ", ".join(absent) + "); "
            "run ./install.sh, whose data step unpacks them"
        )
    work.mkdir(parents=True, exist_ok=True)
    # Written afresh every time: the simulator skips its work when its marker
    # says the reads are there, and the plain files it means are gone (below).
    (work / "selftest.done").unlink(missing_ok=True)
    plain = write_community(specs, out_prefix=work / "selftest", n_pairs=n_pairs, seed=seed)[:2]
    zipped: list[Path] = []
    for path in plain:
        target = path.with_name(path.name + ".gz")
        with path.open("rb") as src, gzip.open(target, "wb", compresslevel=1) as dst:
            shutil.copyfileobj(src, dst, 1 << 20)
        path.unlink()
        zipped.append(target)
    return zipped[0], zipped[1]


def run_lane(lane: str, *, r1: Path, r2: Path, refs_dir: Path, work: Path, threads: int) -> LaneOutcome:
    """One lane, through the code path `openbiota run` uses for it."""
    label = LANES[lane]
    quiet = Reporter(verbose=False)
    t0 = time.monotonic()

    def judged(lineages: Iterable[str]) -> LaneOutcome:
        found = genera_found(lineages)
        return LaneOutcome(lane, label, bool(found), time.monotonic() - t0, "" if found else NOT_FOUND, found)

    try:
        if lane == "genome_profile":
            from openbiota.engines import sylph as sylph_mod  # noqa: PLC0415

            db = refs_dir / "sylph"
            if not sylph_mod.available(db):
                return LaneOutcome(lane, label, False, 0.0, "sylph or its GTDB sketch is not installed - run ./install.sh")
            profile = sylph_mod.run_sylph(sample="selftest", r1=r1, r2=r2, db_dir=db,
                                          work_dir=work / "genome", threads=threads)
            return judged(hit.lineage for hit in (profile.confident_hits if profile is not None else ()))
        from openbiota.expansion import lanes as lanes_mod  # noqa: PLC0415

        block = lanes_mod.run_all(sample="selftest", r1=r1, r2=r2, refs_dir=refs_dir, out_dir=work,
                                  threads=threads, reporter=quiet, only=(lane,))
    except Exception as exc:  # noqa: BLE001 - a self-test reports, it does not crash
        return LaneOutcome(lane, label, False, time.monotonic() - t0, plain_reason(f"{type(exc).__name__}: {exc}"))
    status: dict[str, Any] = (block.get("lane_status") or {}).get(lane) or {}
    state = status.get("state")
    if state == "ran":
        observations = ((block.get("lanes") or {}).get(lane) or {}).get("observations") or []
        return judged(str(o.get("lineage", "")) for o in observations)
    if state == "not_assessed":
        return LaneOutcome(lane, label, False, time.monotonic() - t0,
                           "its tool or reference is not installed - run ./install.sh")
    detail = f"{status.get('reason', '')}\n{status.get('error_detail', '')}"
    return LaneOutcome(lane, label, False, time.monotonic() - t0, plain_reason(detail))


def _measured(lane: str, kwargs: dict[str, Any]) -> LaneOutcome:
    """Runs in a process of its own, so the peak it reads is this lane's and no other's."""
    import resource  # noqa: PLC0415

    from openbiota import projectpath  # noqa: PLC0415

    projectpath.ensure()
    outcome = run_lane(lane, **kwargs)
    # The largest program the lane started: bytes on macOS, kilobytes on Linux.
    peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss * (1 if sys.platform == "darwin" else 1024)
    return replace(outcome, peak_gb=round(peak / 1e9, 1))


def measure_lane(lane: str, **kwargs: Any) -> LaneOutcome:
    """`run_lane`, with the peak memory of the programs it started.

    The operating system reports the largest child a process has ever had,
    not the largest since some moment, so each lane gets a fresh process to
    be the parent of its tools.
    """
    import multiprocessing  # noqa: PLC0415
    from concurrent.futures import ProcessPoolExecutor  # noqa: PLC0415
    from concurrent.futures.process import BrokenProcessPool  # noqa: PLC0415

    t0 = time.monotonic()
    try:
        with ProcessPoolExecutor(max_workers=1, mp_context=multiprocessing.get_context("spawn")) as pool:
            return pool.submit(_measured, lane, kwargs).result()
    except BrokenProcessPool:
        from openbiota.logging_util import OUT_OF_MEMORY  # noqa: PLC0415

        return LaneOutcome(lane, LANES[lane], False, time.monotonic() - t0, OUT_OF_MEMORY)


def run(*, lanes: Sequence[str] | None = None, refs_dir: Path | None = None, threads: int = 8,
        keep: Path | None = None, out: Any = None) -> int:
    """Run the self-test; 0 when every requested lane works."""
    out = out or sys.stdout
    terminal = bool(getattr(out, "isatty", lambda: False)())
    style = Style(terminal and not os.environ.get("NO_COLOR"))
    refs = refs_dir or (REPO / "refs")
    wanted = list(lanes) if lanes else list(LANES)
    unknown = [lane for lane in wanted if lane not in LANES]
    if unknown:
        print(f"unknown lane(s): {', '.join(unknown)}; known: {', '.join(LANES)}", file=sys.stderr)
        return 2

    work = keep or Path(tempfile.mkdtemp(prefix="openbiota-selftest-"))
    work.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    try:
        print(f"Self-test: {len(wanted)} detection lane{'s' if len(wanted) != 1 else ''} on {N_PAIRS:,} made-up "
              f"read pairs from {len(GENOMES)} gut genomes that ship with the repository.\n"
              "A lane passes when it runs to the end and finds the organisms the reads were drawn from.",
              file=out, flush=True)
        try:
            r1, r2 = simulate_reads(work)
        except (OSError, ValueError) as exc:
            print(f"  {style.red('✗')} could not make the test reads: {exc}", file=out)
            return 1
        width = max(len(LANES[lane]) for lane in wanted)
        outcomes: list[LaneOutcome] = []
        redraw = terminal and os.environ.get("TERM", "") != "dumb"
        for lane in wanted:
            # Say which lane is running: the longest takes six minutes. On a
            # terminal the line is replaced by its result; into a log (the
            # installer's, whose status line shows the log's last line) it stays.
            print(f"  · {LANES[lane]} …", end="" if redraw else "\n", file=out, flush=True)
            outcome = measure_lane(lane, r1=r1, r2=r2, refs_dir=refs, work=work, threads=threads)
            if redraw:
                print("\r\033[2K", end="", file=out)
            outcomes.append(outcome)
            took = human_duration(outcome.seconds)
            if outcome.ok:
                memory = f"  {outcome.peak_gb:>4.0f} GB of memory" if outcome.peak_gb else ""
                print(f"  {style.green('✓')} {outcome.label:<{width}}  {took:>7}{memory}", file=out, flush=True)
            else:
                print(f"  {style.red('✗')} {outcome.label:<{width}}  {outcome.reason}", file=out, flush=True)
        failed = [o for o in outcomes if not o.ok]
        total = human_duration(time.monotonic() - started)
        if failed:
            one = len(failed) == 1
            print(f"\n{len(failed)} of {len(outcomes)} lane{'s' if len(outcomes) != 1 else ''} "
                  f"{'does' if one else 'do'} not work on this machine ({total}). "
                  f"A report made now would be missing {'it' if one else 'them'}.", file=out)
            return 1
        print(f"\nEvery lane works ({total}).", file=out)
        return 0
    finally:
        if keep is None:
            shutil.rmtree(work, ignore_errors=True)
