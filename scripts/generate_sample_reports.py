#!/usr/bin/env python3
"""Regenerate the two published sample reports and the documentation screenshots.

    scripts/generate_sample_reports.py [--threads N] [--lp-workers N] [--only ID] [--render-only]

Two real samples ship with the project so a reader can see what a report looks
like before sequencing anything: one community the index places as broadly
disturbed, one it places in the healthy range. Run this whenever the report
changes (a new version, a new section, a changed scale) so the published PDFs
and the pictures in the documentation show the current software.

For each sample the script runs the pipeline on its reads, copies the PDF to
sample-reports/, and renders the pages the documentation shows (summary,
organisms that need attention, what you can do about it, estimated age, disease
patterns) to docs/images/<short id>/. Pages are found by their titles in the
PDF outline, not by page number, so the pictures follow the report when its
layout moves.

Samples run one after another: the pipeline saturates the machine on its own.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OPENBIOTA = REPO / ".venv" / "bin" / "openbiota"
FASTQ_DIR = REPO / "fastq"
RESULTS_DIR = REPO / "results"
OUT_DIR = REPO / "sample-reports"
IMAGES_DIR = REPO / "docs" / "images"

# (sample id, short name used for the image folder, extra flags for the run)
SAMPLES: tuple[tuple[str, str, tuple[str, ...]], ...] = (
    ("EM1_AMD614", "em1", ("--strict-cdiff",)),
    ("MM1_FXX745", "mm1", ()),
)

# Outline title the page is found by (case-insensitive prefix) -> image file name.
# The first level-1 outline entry whose title starts with the key is rendered.
PAGES: tuple[tuple[str, str], ...] = (
    ("Summary", "summary-page.png"),
    ("Organisms that need attention", "organisms-page.png"),
    ("What you can do about it", "evidence-page.png"),
    ("Estimated biological age", "age-page.png"),
    ("Resemblance to published disease patterns at a glance", "profile-page.png"),
)
RENDER_DPI = 144  # A4 at 144 dpi is 1191 x 1684 px: readable in a README, small enough to keep in git


def run_sample(sample: str, extra: tuple[str, ...], *, threads: int, lp_workers: int) -> Path:
    r1, r2 = FASTQ_DIR / f"{sample}_1.fastq.gz", FASTQ_DIR / f"{sample}_2.fastq.gz"
    if not r1.is_file():
        sys.exit(f"{r1} not found: the sample reports are generated from the project's own reads")
    cmd = [str(OPENBIOTA), "run", "--r1", str(r1), "--sample", sample,
           "--threads", str(threads), "--lp-workers", str(lp_workers), *extra]
    if r2.is_file():
        cmd[4:4] = ["--r2", str(r2)]
    log = RESULTS_DIR / sample / "sample_report_run.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    print(f"[{sample}] running the pipeline (log: {log.relative_to(REPO)})", flush=True)
    t0 = time.monotonic()
    with log.open("w") as handle:
        proc = subprocess.run(cmd, cwd=REPO, stdout=handle, stderr=subprocess.STDOUT, check=False)
    minutes = (time.monotonic() - t0) / 60
    if proc.returncode != 0:
        tail = "".join(log.read_text(errors="replace").splitlines(keepends=True)[-25:])
        sys.exit(f"[{sample}] pipeline exited {proc.returncode} after {minutes:.1f} min:\n{tail}")
    pdf = RESULTS_DIR / sample / f"{sample}_report.pdf"
    if not pdf.is_file():
        sys.exit(f"[{sample}] run finished but {pdf.relative_to(REPO)} was not written")
    print(f"[{sample}] done in {minutes:.1f} min", flush=True)
    return pdf


def publish(pdf: Path, sample: str) -> Path:
    OUT_DIR.mkdir(exist_ok=True)
    target = OUT_DIR / f"{sample}_report.pdf"
    shutil.copyfile(pdf, target)
    print(f"[{sample}] {target.relative_to(REPO)} ({target.stat().st_size / 1e6:.1f} MB)")
    return target


def render_pages(pdf: Path, short: str) -> list[Path]:
    try:
        import pymupdf  # the dev extra: pip install -e '.[dev]'
    except ImportError:
        sys.exit("rendering needs PyMuPDF: .venv/bin/pip install -e '.[dev]'")
    doc = pymupdf.open(pdf)
    outline = [(title, page) for level, title, page in doc.get_toc() if level == 1]
    folder = IMAGES_DIR / short
    folder.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for key, filename in PAGES:
        # outline titles carry their section number: "5  Organisms that need attention: ..."
        page_no = next((page for title, page in outline
                        if title.lstrip("0123456789. ").lower().startswith(key.lower())), None)
        if page_no is None:
            print(f"  ! no outline entry starts with {key!r}; {filename} not rendered", file=sys.stderr)
            continue
        pix = doc[page_no - 1].get_pixmap(dpi=RENDER_DPI, alpha=False)
        target = folder / filename
        pix.save(target)
        written.append(target)
        print(f"  {target.relative_to(REPO)}  (page {page_no}, {pix.width}x{pix.height})")
    doc.close()
    return written


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--lp-workers", type=int, default=14)
    ap.add_argument("--only", choices=[s for s, _, _ in SAMPLES], help="just this sample")
    ap.add_argument("--render-only", action="store_true",
                    help="skip the pipeline; publish and render from the PDFs already in results/")
    args = ap.parse_args()

    for sample, short, extra in SAMPLES:
        if args.only and sample != args.only:
            continue
        if args.render_only:
            pdf = RESULTS_DIR / sample / f"{sample}_report.pdf"
            if not pdf.is_file():
                sys.exit(f"{pdf.relative_to(REPO)} does not exist; run without --render-only")
        else:
            pdf = run_sample(sample, extra, threads=args.threads, lp_workers=args.lp_workers)
        publish(pdf, sample)
        render_pages(pdf, short)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
