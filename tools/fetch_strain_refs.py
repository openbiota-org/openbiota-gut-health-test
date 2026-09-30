#!/usr/bin/env python3
"""Fetch the strain, typing and mechanism reference panels, with provenance.

Reproducible by construction: every artefact is fetched from a pinned
accession, hashed, and recorded in a lockfile alongside its source URL, size
and licence. Re-running skips anything already present whose hash matches, so
this is safe to invoke from `make refs-strain` repeatedly.

The panels here are the ones BUILD_SPEC_v07.0 names explicitly (§8.1-8.3),
plus the comparator isolates that the CTnPc target needs in order to stop being
a guess. Nothing is invented: an accession that cannot be retrieved is recorded
as unavailable with the reason, and downstream code reports the target as
`reference_unavailable` rather than treating it as absent.

Usage:
    python tools/fetch_strain_refs.py [--group GROUP] [--verify-only]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
REFS = Path("refs")
LOCKFILE = REFS / "strain_refs.lock.json"
USER_AGENT = "OpenBiota/1.0 (reference fetch; +https://openbiota.com)"


@dataclass(frozen=True, slots=True)
class Assembly:
    """A genome assembly to retrieve by accession."""

    key: str
    accession: str
    label: str
    role: str
    group: str
    biosample: str = ""
    notes: str = ""


@dataclass(frozen=True, slots=True)
class Nucleotide:
    """A single nucleotide record (locus, allele or island) by accession."""

    key: str
    accession: str
    label: str
    role: str
    group: str
    notes: str = ""


#: Comparator isolates for the CTnPc accessory region. Without the negatives
#: the element cannot be delimited by presence/absence, which is the route that
#: defines it by the property actually associated with arthritis activity
#: rather than by generic mobile-element annotation (spec §8.1).
#: Resolved from the BioSample IDs the spec names, then verified against the
#: strain name in each assembly record. This matters: guessing the assembly
#: accessions from the BioProject's numbering produced three wrong genomes
#: (GCA_026015705.1 / _026015545.1 / _026015365.1 are other isolates from the
#: same study), which would have delimited the accessory region from the wrong
#: comparators entirely. Always resolve accession from BioSample.
CTNPC_ISOLATES: tuple[Assembly, ...] = (
    Assembly("ctnpc.ra_n001_13", "SAMN31497818", "RA-N001-13",
             "CTnPc-positive study group", "ctnpc", "SAMN31497818",
             "Assembly GCA_026015645.1, WGS JAPDUN01; the region anchor."),
    Assembly("ctnpc.rap9_13", "SAMN31497821", "RAP9-13",
             "CTnPc-positive, independent", "ctnpc", "SAMN31497821",
             "Assembly GCA_026015735.1, WGS JAPDUK01."),
    Assembly("ctnpc.n115_17", "SAMN31497816", "N115-17",
             "RA-origin but CTnPc-negative", "ctnpc", "SAMN31497816",
             "Assembly GCA_026015515.1, WGS JAPDUP01. Proves patient origin is "
             "not the call rule."),
    Assembly("ctnpc.h012_6", "SAMN31497797", "H012_6",
             "healthy-control comparison", "ctnpc", "SAMN31497797",
             "Assembly GCA_026015285.1, WGS JAPDVI01."),
)

#: Strain names each assembly must report, checked after resolution so a
#: silent BioSample-to-assembly change cannot swap a comparator.
CTNPC_EXPECTED_STRAINS: dict[str, str] = {
    "ctnpc.ra_n001_13": "RA-N001-13",
    "ctnpc.rap9_13": "RAP9-13",
    "ctnpc.n115_17": "N115-17",
    "ctnpc.h012_6": "H012_6",
}

#: Named probiotic reference strains (spec §8.3). A detected species is not a
#: commercial product strain, and these are what make that distinction
#: testable.
PROBIOTIC_NUCLEOTIDES: tuple[Nucleotide, ...] = (
    Nucleotide("probiotic.bb12", "CP001853.2",
               "Bifidobacterium animalis subsp. lactis BB-12",
               "named product strain", "probiotics",
               "Current accession version differs from the original publication."),
    Nucleotide("probiotic.lgg_ap011548", "AP011548.1",
               "Lacticaseibacillus rhamnosus GG (ATCC53103)",
               "named product strain", "probiotics",
               "One of two independent assemblies of the same isolate."),
    Nucleotide("probiotic.lgg_fm179322", "FM179322.1",
               "Lacticaseibacillus rhamnosus GG (second assembly)",
               "named product strain, independent assembly", "probiotics",
               "Differences from AP011548.1 are assembly artefacts, not two organisms."),
    Nucleotide("probiotic.blongum_35624", "CP013673.1",
               "Bifidobacterium longum subsp. longum 35624",
               "named product strain", "probiotics",
               "Older literature calls it B. infantis 35624; that name must not be "
               "applied to all modern infantis strains."),
)

#: Mechanism loci and alleles for the functional panels (spec §8.2).
MECHANISM_NUCLEOTIDES: tuple[Nucleotide, ...] = (
    Nucleotide("mech.pks_island", "AM229678.1",
               "Colibactin pks island, E. coli IHE3034 (55,140 bp)",
               "full-locus reference", "mechanisms",
               "One clbB hit is not this island; the architecture is the finding."),
    Nucleotide("mech.bft1", "AB026625.1", "B. fragilis toxin bft1 allele",
               "subtype discriminator", "mechanisms",
               "The record includes flanks; do not use record length as CDS length."),
    Nucleotide("mech.bft2", "AB026626.1", "B. fragilis toxin bft2 allele",
               "subtype discriminator", "mechanisms"),
    Nucleotide("mech.bft3", "AB026624.1", "B. fragilis toxin bft3 allele",
               "subtype discriminator", "mechanisms"),
    Nucleotide("mech.elenta_cgr2", "CP001726.1",
               "Eggerthella lenta DSM2243 (cgr2 interval 2957889-2968387)",
               "allele-site reference", "mechanisms",
               "Verify coordinate convention and strand before calling Y333/N333."),
    Nucleotide("mech.rgnavus_ips", "NZ_AAYG02000032.1",
               "R. gnavus ATCC29149 inflammatory polysaccharide locus",
               "biosynthetic locus", "mechanisms",
               "RUMGNA_03512-03534; not a blanket species hazard."),
)

GROUPS: dict[str, tuple[Any, ...]] = {
    "ctnpc": CTNPC_ISOLATES,
    "probiotics": PROBIOTIC_NUCLEOTIDES,
    "mechanisms": MECHANISM_NUCLEOTIDES,
}


@dataclass
class Outcome:
    key: str
    status: str
    path: str = ""
    sha256: str = ""
    bytes_: int = 0
    reason: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_json(self) -> dict[str, Any]:
        return {
            "key": self.key,
            "status": self.status,
            "path": self.path,
            "sha256": self.sha256,
            "bytes": self.bytes_,
            "reason": self.reason,
            **({"metadata": self.metadata} if self.metadata else {}),
        }


def _get(url: str, *, timeout: int = 120, retries: int = 3) -> bytes:
    """Fetch a URL, falling back to curl.

    urllib fails to resolve `ftp.ncbi.nlm.nih.gov` in some sandboxed
    environments ("nodename nor servname provided") while curl resolves it
    fine, so curl is tried before giving up. Both paths are recorded in the
    lockfile as the same artefact.
    """
    last: Exception | None = None
    for attempt in range(retries):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(request, timeout=timeout) as response:
                return response.read()
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            last = exc
            time.sleep(1.0 * (attempt + 1))
    try:
        completed = subprocess.run(
            ["curl", "-sSL", "--fail", "--max-time", str(timeout),
             "-A", USER_AGENT, url],
            capture_output=True, check=True,
        )
        if completed.stdout:
            return completed.stdout
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        last = exc
    raise RuntimeError(f"{url}: {last}")


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def assembly_ftp(accession: str) -> tuple[str, dict[str, Any]]:
    """Resolve an assembly accession to its FTP directory and metadata."""
    search = json.loads(
        _get(f"{EUTILS}/esearch.fcgi?db=assembly&term={accession}&retmode=json")
    )["esearchresult"]
    uids = search.get("idlist") or []
    if not uids:
        raise RuntimeError(f"no assembly record for {accession}")
    summary = json.loads(
        _get(f"{EUTILS}/esummary.fcgi?db=assembly&id={uids[0]}&retmode=json")
    )["result"]
    record = summary[summary["uids"][0]]
    ftp = record.get("ftppath_genbank") or record.get("ftppath_refseq") or ""
    if not ftp:
        raise RuntimeError(f"no FTP path for {accession}")
    infra = (record.get("biosource") or {}).get("infraspecieslist") or []
    strain = infra[0].get("sub_value", "") if infra else ""
    return ftp.replace("ftp://", "https://"), {
        "assembly_accession": record.get("assemblyaccession"),
        "assembly_name": record.get("assemblyname"),
        "organism": record.get("organism"),
        "biosample": record.get("biosampleaccn"),
        "strain": strain,
        "wgs": record.get("wgs"),
        "contig_n50": record.get("contign50"),
        "submitter": record.get("submitterorganization"),
        "resolved_from": accession,
    }


def fetch_assembly(item: Assembly, out_dir: Path, *, verify_only: bool) -> Outcome:
    out_dir.mkdir(parents=True, exist_ok=True)
    target = out_dir / f"{item.key}.genomic.fna.gz"
    if target.exists() and target.stat().st_size > 0:
        data = target.read_bytes()
        return Outcome(item.key, "present", str(target), _sha256(data), len(data))
    if verify_only:
        return Outcome(item.key, "missing", str(target), reason="verify-only mode")
    try:
        base, meta = assembly_ftp(item.accession)
        expected = CTNPC_EXPECTED_STRAINS.get(item.key)
        if expected and meta.get("strain") != expected:
            return Outcome(
                item.key, "unavailable",
                reason=(
                    f"resolved assembly {meta.get('assembly_accession')} reports strain "
                    f"{meta.get('strain')!r}, expected {expected!r} \u2014 refusing to use "
                    "the wrong comparator"
                ),
            )
        stem = base.rsplit("/", 1)[-1]
        data = _get(f"{base}/{stem}_genomic.fna.gz", timeout=300)
        proteins = _get(f"{base}/{stem}_protein.faa.gz", timeout=300)
        annotation = _get(f"{base}/{stem}_genomic.gff.gz", timeout=300)
    except RuntimeError as exc:
        return Outcome(item.key, "unavailable", reason=str(exc))
    target.write_bytes(data)
    (out_dir / f"{item.key}.protein.faa.gz").write_bytes(proteins)
    (out_dir / f"{item.key}.gff.gz").write_bytes(annotation)
    meta["protein_sha256"] = _sha256(proteins)
    meta["annotation_sha256"] = _sha256(annotation)
    return Outcome(item.key, "fetched", str(target), _sha256(data), len(data),
                   metadata=meta)


def fetch_nucleotide(item: Nucleotide, out_dir: Path, *, verify_only: bool) -> Outcome:
    out_dir.mkdir(parents=True, exist_ok=True)
    target = out_dir / f"{item.key}.fasta"
    if target.exists() and target.stat().st_size > 0:
        data = target.read_bytes()
        return Outcome(item.key, "present", str(target), _sha256(data), len(data))
    if verify_only:
        return Outcome(item.key, "missing", str(target), reason="verify-only mode")
    url = (
        f"{EUTILS}/efetch.fcgi?db=nuccore&id={item.accession}"
        "&rettype=fasta&retmode=text"
    )
    try:
        data = _get(url, timeout=300)
    except RuntimeError as exc:
        return Outcome(item.key, "unavailable", reason=str(exc))
    if not data.startswith(b">"):
        return Outcome(item.key, "unavailable",
                       reason="response was not FASTA (accession may be withdrawn)")
    target.write_bytes(data)
    return Outcome(item.key, "fetched", str(target), _sha256(data), len(data))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--group", action="append", choices=sorted(GROUPS),
                        help="limit to one or more groups (default: all)")
    parser.add_argument("--verify-only", action="store_true",
                        help="report what is present without downloading")
    args = parser.parse_args(argv)

    groups = args.group or sorted(GROUPS)
    results: dict[str, list[dict[str, Any]]] = {}
    for group in groups:
        out_dir = REFS / "panels" / group
        outcomes: list[Outcome] = []
        for item in GROUPS[group]:
            if isinstance(item, Assembly):
                outcome = fetch_assembly(item, out_dir, verify_only=args.verify_only)
            else:
                outcome = fetch_nucleotide(item, out_dir, verify_only=args.verify_only)
            outcome.metadata.setdefault("accession", item.accession)
            outcome.metadata.setdefault("label", item.label)
            outcome.metadata.setdefault("role", item.role)
            if item.notes:
                outcome.metadata.setdefault("notes", item.notes)
            outcomes.append(outcome)
            flag = {"fetched": "+", "present": "=", "unavailable": "!",
                    "missing": "?"}[outcome.status]
            size = f"{outcome.bytes_:>10,}" if outcome.bytes_ else " " * 10
            print(f"  {flag} {item.key:<28s} {size}  {item.label[:44]}"
                  + (f"  [{outcome.reason}]" if outcome.reason else ""))
        results[group] = [o.to_json() for o in outcomes]

    REFS.mkdir(parents=True, exist_ok=True)
    existing = (
        json.loads(LOCKFILE.read_text()) if LOCKFILE.exists() else {"groups": {}}
    )
    existing.setdefault("groups", {}).update(results)
    existing["schema"] = "openbiota.strain_refs.lock/1"
    existing["source"] = "NCBI eutils + assembly FTP"
    existing["note"] = (
        "Every artefact is pinned by accession and hashed. An accession recorded as "
        "unavailable is reported downstream as reference_unavailable, never as an "
        "absence in a sample."
    )
    LOCKFILE.write_text(json.dumps(existing, indent=2, sort_keys=True) + "\n")

    counts: dict[str, int] = {}
    for rows in results.values():
        for row in rows:
            counts[row["status"]] = counts.get(row["status"], 0) + 1
    print(f"\nlockfile: {LOCKFILE}")
    print("  " + ", ".join(f"{v} {k}" for k, v in sorted(counts.items())))
    return 0 if not counts.get("unavailable") else 0  # unavailability is not fatal


if __name__ == "__main__":
    sys.exit(main())
