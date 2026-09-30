"""Release-aware crosswalks from MetaPhlAn SGBs to GTDB R232 species.

Why this is built from assemblies and not from names
----------------------------------------------------
MetaPhlAn ships a bridge from each of its databases to *one* GTDB release:
Jun23 -> R207, Jan26 -> R226. There is no R232 bridge, and this repository
had been applying the Jan25/R220 file to Jun23 profiles and reading the
result beside an R232 Sylph lane - two release mismatches at once. A GTDB
species name is not stable across releases: species are renamed, split
and merged as genomes are added. So the bridge's R226 name is used only as
the first step, and the walk onward to R232 goes through the genomes:

    SGB  --bridge-->  R226 species  --sp_clusters_r226-->  member accessions
         --R232 metadata-->  the R232 species those same genomes now sit in

What the members say decides the relationship, and the relationship is
stored, not flattened:

    exact       every member still in one R232 species with the same name
    synonym     every member in one R232 species under a different name
    split       members now in several R232 species  -> keep as a complex
    merge       several SGBs' members land in one R232 species
    parent      bridge gives no species; genus (or higher) is what is known
    unresolved  the bridge names a species R226 no longer lists, or all
                members are absent from R232

A `split` is rendered as one unresolved complex, never as several certain
species. Two SGBs in a `merge` are two native observations of one canonical
organism and render once. Numeric SGB suffixes are never compared across
releases.

Outputs: `refs/crosswalk/<db>_to_gtdb_r232.tsv.gz` plus a JSON summary
with counts per relationship, and a fingerprint of the inputs so a change
to any of them invalidates the file.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Final

from openbiota.expansion import gtdb

CROSSWALK_DIR: Final = Path("refs/crosswalk")
UPSTREAM_DIR: Final = CROSSWALK_DIR / "upstream"

#: database -> (upstream bridge file, GTDB release the bridge targets)
BRIDGES: Final[dict[str, tuple[str, str]]] = {
    "mpa_vJan26_CHOCOPhlAnSGB_202605": ("mpa_vJan26_CHOCOPhlAnSGB_202605_SGB2GTDB_r226.tsv", "r226"),
    "mpa_vJun23_CHOCOPhlAnSGB_202403": ("mpa_vJun23_CHOCOPhlAnSGB_202403_SGB2GTDB_r207.tsv", "r207"),
}
TARGET: Final = "r232"
RELATIONSHIPS: Final = ("exact", "synonym", "split", "merge", "parent", "unresolved")


@dataclass(frozen=True)
class Edge:
    sgb: str
    bridge_release: str
    bridge_species: str          # R226/R207 species the bridge names ('' when none)
    bridge_taxonomy: str
    relationship: str
    canonical_species: str       # R232 species when exact/synonym; '' otherwise
    canonical_genus: str         # best genus known in R232 terms (or bridge genus)
    r232_species_members: str    # "s__A:12;s__B:3" - the evidence
    n_members_checked: int
    n_members_in_r232: int
    complex_members: str         # for split: the R232 species in the complex, ';'-joined


def _bridge(db: str) -> tuple[dict[str, str], str]:
    fname, release = BRIDGES[db]
    path = UPSTREAM_DIR / fname
    out: dict[str, str] = {}
    with path.open() as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[0].startswith("SGB"):
                out[parts[0]] = parts[1]
    return out, release


def _fingerprint(paths: list[Path]) -> str:
    h = hashlib.sha256()
    for p in paths:
        h.update(p.name.encode())
        h.update(str(p.stat().st_size).encode())
    return h.hexdigest()[:16]


def build(db: str, *, force: bool = False) -> Path:
    """Build (or reuse) the crosswalk for one MetaPhlAn database."""
    fname, bridge_release = BRIDGES[db]
    out = CROSSWALK_DIR / f"{db}_to_gtdb_{TARGET}.tsv.gz"
    summary_path = out.with_suffix("").with_suffix(".summary.json")
    inputs = [UPSTREAM_DIR / fname, gtdb.GTDB_DIR / bridge_release / f"sp_clusters_{bridge_release}.tsv",
              gtdb.build_slim(TARGET)]
    fp = _fingerprint(inputs)
    if out.is_file() and summary_path.is_file() and not force \
            and json.loads(summary_path.read_text()).get("input_fingerprint") == fp:
        return out

    bridge, _ = _bridge(db)
    clusters = gtdb.species_clusters(bridge_release)
    r232 = gtdb.genome_table(TARGET)

    edges: list[Edge] = []
    # First pass: per-SGB resolution through members.
    by_canonical: dict[str, list[str]] = defaultdict(list)
    for sgb, taxonomy in bridge.items():
        sp = gtdb.species_of(taxonomy)
        genus = gtdb.genus_of(taxonomy)
        if not sp:
            edges.append(Edge(sgb, bridge_release, "", taxonomy, "parent", "", genus, "", 0, 0, ""))
            continue
        cluster = clusters.get(sp)
        if cluster is None:
            edges.append(Edge(sgb, bridge_release, sp, taxonomy, "unresolved", "", genus, "", 0, 0, ""))
            continue
        # Where did this species' genomes go in R232?
        counts: Counter[str] = Counter()
        seen = 0
        for acc in cluster.members:
            g = r232.get(acc)
            if g is None:
                # try unversioned match (a version bump keeps the genome)
                continue
            seen += 1
            if g.species:
                counts[g.species] += 1
        checked = len(cluster.members)
        if not counts:
            edges.append(Edge(sgb, bridge_release, sp, taxonomy, "unresolved", "", genus, "", checked, seen, ""))
            continue
        evidence = ";".join(f"{k}:{v}" for k, v in counts.most_common())
        top, top_n = counts.most_common(1)[0]
        total = sum(counts.values())
        # A split is real when a second R232 species holds a material share
        # of the members - not one stray reassigned genome in a thousand.
        others = [(k, v) for k, v in counts.items() if k != top]
        material = [k for k, v in others if v >= max(2, 0.05 * total)]
        if material:
            complex_members = ";".join([top, *material])
            edges.append(Edge(sgb, bridge_release, sp, taxonomy, "split", "", gtdb.genus_of(
                next((g.taxonomy for g in (r232.get(a) for a in cluster.members) if g and g.species == top), taxonomy)),
                evidence, checked, seen, complex_members))
            continue
        canonical_genus = ""
        for acc in cluster.members:
            g = r232.get(acc)
            if g and g.species == top:
                canonical_genus = gtdb.genus_of(g.taxonomy)
                break
        rel = "exact" if top == sp else "synonym"
        edges.append(Edge(sgb, bridge_release, sp, taxonomy, rel, top, canonical_genus or genus, evidence, checked, seen, ""))
        by_canonical[top].append(sgb)

    # Second pass: merges - several SGBs onto one R232 species.
    merged_into: dict[str, str] = {}
    for canonical, sgbs in by_canonical.items():
        if len(sgbs) > 1:
            for s in sgbs:
                merged_into[s] = canonical
    final: list[Edge] = []
    for e in edges:
        if e.sgb in merged_into and e.relationship in ("exact", "synonym"):
            final.append(Edge(**{**asdict(e), "relationship": "merge"}))
        else:
            final.append(e)

    CROSSWALK_DIR.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".tmp")
    with gzip.open(tmp, "wt", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([f.name for f in Edge.__dataclass_fields__.values()])
        for e in sorted(final, key=lambda e: (int(e.sgb[3:]) if e.sgb[3:].isdigit() else 0)):
            w.writerow([getattr(e, f) for f in Edge.__dataclass_fields__])
    tmp.replace(out)
    counts_by_rel = Counter(e.relationship for e in final)
    summary = {
        "database": db, "bridge_release": bridge_release, "target_release": TARGET,
        "n_sgbs": len(final), "relationships": dict(counts_by_rel),
        "n_canonical_species": len({e.canonical_species for e in final if e.canonical_species}),
        "n_merge_groups": sum(1 for v in by_canonical.values() if len(v) > 1),
        "input_fingerprint": fp,
        "inputs": [str(p) for p in inputs],
        "method": (
            "bridge species -> member accessions of that species in the bridge release's sp_clusters -> "
            "the R232 species those accessions now carry; relationship typed from the distribution; a second "
            "R232 species holding >=5% (and >=2) of members makes a split; several SGBs on one R232 species "
            "make a merge"
        ),
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    return out


def load(db: str) -> dict[str, Edge]:
    """sgb -> Edge for a built crosswalk (builds it if needed)."""
    path = build(db)
    out: dict[str, Edge] = {}
    with gzip.open(path, "rt") as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            row["n_members_checked"] = int(row["n_members_checked"])
            row["n_members_in_r232"] = int(row["n_members_in_r232"])
            e = Edge(**row)
            out[e.sgb] = e
    return out


def summary(db: str) -> dict:
    build(db)
    return json.loads((CROSSWALK_DIR / f"{db}_to_gtdb_{TARGET}.summary.json").read_text())


if __name__ == "__main__":  # pragma: no cover
    for name in BRIDGES:
        build(name)
        print(json.dumps(summary(name), indent=1))


def release_species_map(from_release: str, *, force: bool = False) -> dict[str, str]:
    """Every species name in an older GTDB release -> its R232 species.

    Walked by genomes, as the SGB crosswalk is: the old species' member
    accessions looked up in R232, majority species taken when it holds at
    least half the members that survive. A species whose members scatter or
    vanish maps to '' and stays in its own release's name, flagged.
    Cached to refs/crosswalk/gtdb_<release>_to_r232.species.tsv.gz.
    """
    out_path = CROSSWALK_DIR / f"gtdb_{from_release}_to_{TARGET}.species.tsv.gz"
    if out_path.is_file() and not force:
        result: dict[str, str] = {}
        with gzip.open(out_path, "rt") as fh:
            for line in fh:
                a, _, b = line.rstrip("\n").partition("\t")
                if a:
                    result[a] = b
        return result
    clusters = gtdb.species_clusters(from_release)
    r232 = gtdb.genome_table(TARGET)
    by_unv = {gtdb.unversioned(a): g for a, g in r232.items()}
    result = {}
    for sp, cluster in clusters.items():
        counts: Counter[str] = Counter()
        for acc in cluster.members[:2000]:
            g = r232.get(acc) or by_unv.get(gtdb.unversioned(acc))
            if g is not None and g.species:
                counts[g.species] += 1
        if counts:
            top, n = counts.most_common(1)[0]
            result[sp] = top if n >= 0.5 * sum(counts.values()) else ""
        else:
            result[sp] = ""
    CROSSWALK_DIR.mkdir(parents=True, exist_ok=True)
    with gzip.open(out_path, "wt") as fh:
        for sp, dest in sorted(result.items()):
            fh.write(f"{sp}\t{dest}\n")
    return result
