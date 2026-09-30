#!/usr/bin/env python3
"""Write docs/ATTRIBUTION.md: every tool and reference the report is built on, with its licence.

The project redistributes none of these; each is fetched from its publisher
at install time and pinned by version and checksum under `refs/expanded/`.
Attribution is nonetheless owed, and a reader of the report is entitled to
know what the numbers rest on. Reference licences are read from the lock
files; tool pins from `refs/expanded/tools.lock.json` where it exists (it
is written per machine and not tracked), else from the pinned versions
below.
"""

from __future__ import annotations

import datetime as dt
import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
OUT = REPO / "docs" / "ATTRIBUTION.md"
LOCKS = REPO / "refs" / "expanded" / "locks"

#: name -> (what it does here, licence as the project publishes it, home)
TOOLS: dict[str, tuple[str, str, str]] = {
    "DIAMOND": ("translated search of every read against the metabolite gene panels", "GPL-3.0",
                "https://github.com/bbuchfink/diamond"),
    "Bowtie2": ("host-read removal; competitive confirmation alignments", "GPL-3.0",
                "https://github.com/BenLangmead/bowtie2"),
    "samtools": ("alignment sorting, indexing and depth", "MIT/Expat", "https://github.com/samtools/samtools"),
    "MetaPhlAn 3.1": ("the scoring lane the reference cohort and GMWI2 were built with", "MIT",
                      "https://github.com/biobakery/MetaPhlAn"),
    "MetaPhlAn 4 (Jun23, Jan26), StrainPhlAn 4, sample2markers": (
        "species-level genome bins; comparative strain placement", "MIT (software); CHOCOPhlAn databases CC BY-NC-SA",
        "https://github.com/biobakery/MetaPhlAn"),
    "PhyloPhlAn, RAxML": ("marker alignment and tree inference inside StrainPhlAn", "MIT; GPL-3.0",
                          "https://github.com/biobakery/phylophlan"),
    "sylph": ("whole-genome containment against GTDB R232 and GlobDB r232", "MIT",
              "https://github.com/bluenote-1577/sylph"),
    "mOTUs 4.1": ("universal single-copy marker profiling", "GPL-3.0", "https://github.com/motu-tool/mOTUs"),
    "Kraken 2, Bracken": ("read classification against UHGG v2.0.2 and the gut rescue panel", "MIT; GPL-3.0",
                          "https://github.com/DerrickWood/kraken2"),
    "SingleM (with OrfM, smafa, mfqe, HMMER)": ("marker-window profiling of lineages no catalogue names",
                                                "GPL-3.0 (SingleM, HMMER); MIT (OrfM, smafa, mfqe)",
                                                "https://github.com/wwood/singlem"),
    "metaSPAdes": ("targeted assembly of reads no method could place", "GPL-2.0",
                   "https://github.com/ablab/spades"),
    "skani, FastANI": ("genome-wide ANI for supplement reconciliation", "MIT; Apache-2.0",
                       "https://github.com/bluenote-1577/skani"),
    "seqkit": ("read selection for targeted assembly; marker extraction", "MIT",
               "https://github.com/shenwei356/seqkit"),
    "aria2": ("parallel, resumable reference downloads", "GPL-2.0-or-later", "https://aria2.github.io"),
    "fastp": ("read QC and trimming", "MIT", "https://github.com/OpenGene/fastp"),
    "GMWI2": ("the published gut-health index", "MIT", "https://github.com/danielchang2002/GMWI2"),
    "MICOM (with HiGHS, OSQP)": ("model-assisted metabolic scenarios", "Apache-2.0; MIT; Apache-2.0",
                                 "https://github.com/micom-dev/micom"),
    "ReportLab": ("the PDF", "BSD-3-Clause", "https://www.reportlab.com/opensource/"),
    "NumPy, PyYAML, scikit-learn": ("arithmetic, configuration, the age model", "BSD-3-Clause; MIT; BSD-3-Clause",
                                    "https://numpy.org"),
    "curatedMetagenomicData": ("the 3,027-adult reference cohort's profiles and metadata", "Artistic-2.0",
                               "https://waldronlab.io/curatedMetagenomicData/"),
    "AGORA2 / VMH": ("genome-scale metabolic models behind the scenarios", "CC BY 4.0",
                     "https://www.vmh.life"),
    "BacDive (DSMZ)": ("culture-collection identities and reference-strain traits, via API v2",
                       "CC BY 4.0; DSMZ asks to be contacted for commercial use", "https://bacdive.dsmz.de"),
    "UniProt": ("reference proteins for the metabolite gene panels", "CC BY 4.0", "https://www.uniprot.org"),
}


def references() -> list[tuple[str, str, str, str]]:
    rows: list[tuple[str, str, str, str]] = []
    for lock in sorted(LOCKS.glob("*.lock.json")):
        d = json.loads(lock.read_text())
        urls = d.get("urls") or []
        if not urls and isinstance(d.get("files"), dict):
            urls = [v.get("url") for v in d["files"].values() if isinstance(v, dict) and v.get("url")]
        home = urls[0].rsplit("/", 1)[0] if urls else ""
        count = d.get("source_count")
        unit = d.get("count_unit") or ""
        size = f"{count:,} {unit}" if isinstance(count, int) and unit else ""
        rows.append((str(d.get("source_id") or lock.stem.replace(".lock", "")),
                     f"{d.get('release_id') or ''} {('(' + size + ')') if size else ''}".strip(),
                     str(d.get("license") or d.get("licence") or "see source"), home))
    return rows


def main() -> int:
    lines = [
        "# Attribution",
        "",
        f"Generated {dt.date.today().isoformat()} by `scripts/attribution_doc.py`. OpenBiota redistributes none of the",
        "software or data below: each is fetched from its publisher at install time and pinned by version and",
        "checksum under `refs/expanded/`. Every run records the exact versions and release identifiers it used in",
        "`results.json` (`detection.reference_releases`, `run`) and in `run.log`. Licences are as each project publishes",
        "them at the time of writing; check the source before any commercial use, and note that the MetaPhlAn",
        "databases and BacDive carry noncommercial or contact-first terms of their own.",
        "",
        "## Software",
        "",
        "| tool | used for | licence | home |",
        "|---|---|---|---|",
    ]
    for name, (what, lic, home) in TOOLS.items():
        lines.append(f"| {name} | {what} | {lic} | <{home}> |")
    lines += ["", "## Reference releases", "",
              "Each lock file under `refs/expanded/locks/` carries the URLs, sizes, SHA-256 digests and retrieval date.",
              "", "| source | release | licence | home |", "|---|---|---|---|"]
    for sid, rel, lic, home in references():
        lines.append(f"| {sid} | {rel} | {lic} | {('<' + home + '>') if home else ''} |")
    lines += [
        "",
        "## Published research",
        "",
        "Every metabolite pathway, disease-pattern profile, microbial group, organism description and intervention",
        "assertion names the study it rests on (DOI or PubMed identifier) in its definition file and on the page of",
        "the report where it is used. The complete catalogue of sources is on the website's research page and in",
        "`specs/research/`.",
        "",
    ]
    OUT.write_text("\n".join(lines))
    print(f"wrote {OUT} ({len(TOOLS)} tools, {len(references())} reference releases)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
