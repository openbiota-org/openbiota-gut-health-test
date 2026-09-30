"""Reference resolution, acquisition, audit and locking (BUILD_SPEC_v05.0 §3).

Four entry points, matching the spec's required commands:

* `resolve_catalog` — names and aliases to a frozen taxonomy, then to actual
  assembly accessions. Emits every gap.
* `acquire` — download sequence by `accession.version`, checksummed.
* `build_bundle` — mask, compute informative regions competitively, index,
  and write the installed manifest.
* `audit_bundle` — re-derive coverage from the artifacts on disk and fail any
  claim the sequence does not support.

The discipline that matters most here is stated once in the spec and repeated
throughout: **a target is not covered because its name exists in taxonomy.**
Coverage is built from actual accessions, actual unmasked bases and actual
demonstrated ability to distinguish the target from its relatives. So:

* Name resolution is exact. It matches the requested name or a declared alias
  against NCBI's scientific name and synonyms. It never takes the first
  substring match, and an ambiguous name becomes a review record rather than
  a guess (acceptance P019).
* An assembly under `--max-asset-bytes` installs as a genome; one over it
  installs its organelle or marker references instead and records
  `organelle_only`/`marker_only`. That is a real, visible downgrade of what
  the target can claim — never a silent substitution.
* Informative regions are computed *after* everything is loaded together, by
  competing every target against every other target, the host and the food
  and near-neighbour decoys. A region shared with a relative is not
  species-specific evidence, whatever its coverage depth.
* A checksum mismatch invalidates the affected asset and everything derived
  from it, and nothing else (acceptance P020).
"""

from __future__ import annotations

import concurrent.futures as futures
import fcntl
import hashlib
import json
import re
import shutil
import subprocess
import time
import zipfile
from collections.abc import Callable, Mapping, MutableMapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Final

import numpy as np

from openbiota import net
from openbiota.errors import ReferenceFetchError
from openbiota.pathogens.catalog import Catalog, Seed
from openbiota.pathogens.kmers import SCALED, SHARED, KmerIndex, iter_fasta, select_kmers
from openbiota.pathogens.schema import PathogenError, ReferenceLock

__all__ = [
    "BundleManifest",
    "Resolution",
    "acquire",
    "audit_bundle",
    "build_bundle",
    "resolve_catalog",
]

_DATASETS: Final = "https://api.ncbi.nlm.nih.gov/datasets/v2alpha"
_EUTILS: Final = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

#: Canonical genomic bin size for independent-region counting (spec §5.4
#: rule 2). Three adjacent windows of one locus are one region, so regions are
#: fixed bins of the reference rather than arbitrary alignment clusters.
BIN_BASES: Final[int] = 10_000

#: A bin counts as target-specific only when it carries at least this many
#: k-mers seen in no other reference. One incidental unique k-mer in an
#: otherwise conserved bin is not discriminating sequence.
MIN_INFORMATIVE_KMERS_PER_BIN: Final[int] = 3

#: Marker search terms, by group. Organelle and ribosomal loci are the
#: published route for the large-genome eukaryotes; they are installed as
#: marker evidence and can never yield whole-genome breadth (spec §9.7 rule 3).
_MARKER_QUERIES: Final[dict[str, tuple[str, ...]]] = {
    "helminths": (
        "mitochondrion complete genome",
        "internal transcribed spacer",
        "18S ribosomal RNA",
        "cytochrome c oxidase subunit 1",
    ),
    "protozoa": (
        "18S ribosomal RNA",
        "internal transcribed spacer",
        "mitochondrion complete genome",
    ),
    "microsporidia": ("18S ribosomal RNA", "internal transcribed spacer"),
    "fungi": ("internal transcribed spacer", "28S ribosomal RNA"),
    "other_eukaryotes": ("18S ribosomal RNA", "internal transcribed spacer"),
    "bacteria": ("16S ribosomal RNA",),
    "dna_viruses": ("complete genome",),
    "rna_viruses": ("complete genome",),
    "other_viruses": ("complete genome",),
}

#: Assembly levels, best first. A rare target with only a contig assembly is
#: still installed: excluding it would manufacture a reference gap (spec §3.7
#: step 2).
_LEVEL_RANK: Final[dict[str, int]] = {
    "Complete Genome": 0,
    "Chromosome": 1,
    "Scaffold": 2,
    "Contig": 3,
}
_CATEGORY_RANK: Final[dict[str, int]] = {
    "reference genome": 0,
    "representative genome": 1,
}


# --------------------------------------------------------------------------- #
# Resolution
# --------------------------------------------------------------------------- #


@dataclass
class Resolution:
    """What the build could actually find for one seed."""

    seed_id: str
    requested_name: str
    resolution_state: str = "pending_resolution"
    accepted_taxon: str | None = None
    ncbi_taxids: tuple[int, ...] = ()
    accession: str | None = None
    assembly_level: str | None = None
    assembly_bytes: int = 0
    marker_accessions: tuple[str, ...] = ()
    route: str = "none"
    gap_reason: str | None = None
    candidates_considered: tuple[str, ...] = ()
    note: str = ""
    #: The seed whose installed reference already carries this row's sequence.
    covered_by: str | None = None

    def to_json(self) -> dict[str, Any]:
        return {
            "seed_id": self.seed_id,
            "requested_name": self.requested_name,
            "resolution_state": self.resolution_state,
            "accepted_taxon": self.accepted_taxon,
            "ncbi_taxids": list(self.ncbi_taxids),
            "accession": self.accession,
            "assembly_level": self.assembly_level,
            "assembly_bytes": self.assembly_bytes,
            "marker_accessions": list(self.marker_accessions),
            "route": self.route,
            "gap_reason": self.gap_reason,
            "candidates_considered": list(self.candidates_considered),
            "note": self.note,
            "covered_by": self.covered_by,
        }


def _get_json(url: str, *, timeout: float = 90.0) -> Any:
    response = net.get(url, timeout=timeout, attempts=4)
    try:
        return json.loads(response.body)
    except json.JSONDecodeError as exc:
        raise ReferenceFetchError(f"{url}: response was not JSON ({exc})") from exc


def _quote(value: str) -> str:
    from urllib.parse import quote

    return quote(value, safe="")


#: Wording a catalogue row uses to describe a *set* of organisms. NCBI has no
#: taxon called "Fusarium, clinically relevant species complexes", so a row
#: written that way failed to resolve and was reported as never searched. The
#: taxon behind it is the genus, which is the resolution the row already
#: declares, so the qualifier is stripped and the genus searched.
_SET_QUALIFIERS: Final = re.compile(
    r"\b(species complex(es)?|sensu lato|sensu stricto|complex(es)?|group|spp\.?|"
    r"clinically (relevant|documented) (members|species|species complexes)|"
    r"human[- ]infecting (members|species)|verified relevant relatives|relevant relatives|"
    r"recognized genotypes|curated (human[- ]infecting members|dietary species|human subtypes)|"
    r"pathogenic clade|unnamed|subsp\.|subspecies|serovar|serotype|biotype|pathotype|lineage [A-Z])\b",
    re.IGNORECASE,
)


#: Words a catalogue row starts with that are English, not taxonomy. Stripping
#: a row's qualifiers can leave one of these standing alone, and NCBI will
#: happily match it: "Human torovirus-like findings" reduced to "Human" and
#: resolved to *Homo sapiens*, which would have installed the 3.1 Gbp host
#: genome as a virus reference and matched every human read in the sample.
_NOT_TAXON_WORDS: Final[frozenset[str]] = frozenset({
    "human", "humans", "animal", "animals", "nonhuman", "other", "others", "additional",
    "common", "free", "neighboring", "neighbouring", "fish", "reclassified", "curated",
    "clinically", "verified", "recognized", "recognised", "hybrid", "diffusely",
    "enteroaggregative", "enteropathogenic", "enterotoxigenic", "enteroinvasive",
    "glabrata", "food", "plant", "shiga", "supported", "documented", "relevant", "unnamed",
})

#: Host and near-host taxa a pathogen reference must never resolve to. The
#: host genome is already a decoy; installing it as a target would make the
#: target match everything.
_FORBIDDEN_TAXA: Final[frozenset[str]] = frozenset({
    "homo sapiens", "homo", "hominidae", "primates", "mammalia", "chordata", "eukaryota",
    "bacteria", "archaea", "viruses", "fungi", "metazoa", "root", "cellular organisms",
})


def _name_variants(name: str) -> list[str]:
    """Progressively plainer spellings of one catalogue name.

    A catalogue row is written for a reader ("Cyclospora cayetanensis/lineage
    A", "Aspergillus, clinically relevant cryptic species", "Animal Giardia
    competitors"); NCBI wants a taxon. Each variant drops one layer of the
    row's own wording and ends at the organism name embedded in it, which for
    a genus- or group-level row is the reference the row is asking for. Most
    specific first, so a species row never resolves to its genus while the
    species itself is available.
    """
    out: list[str] = []

    def add(value: str) -> None:
        value = re.sub(r"\s+", " ", value).strip(" ,;/-.")
        if value and value not in out:
            out.append(value)

    def taxon_like(word: str) -> bool:
        return (
            len(word) > 3
            and word[:1].isupper()
            and word.isalpha()
            and word.lower() not in _NOT_TAXON_WORDS
        )

    add(name)
    plain = re.sub(r"\([^)]*\)", " ", name)
    add(plain)
    add(re.split(r"[/(]", name)[0])
    head = plain.split(",")[0]
    add(head)
    # A serovar or subspecies row is filed under the subspecies in NCBI, not
    # under the species: "Salmonella enterica serovar Paratyphi A" lives at
    # "Salmonella enterica subsp. enterica serovar Paratyphi A". Tried before
    # anything that would strip the serovar, or the row silently becomes the
    # species and duplicates a reference already installed.
    serovar = re.search(r"\bserovar\s+([A-Za-z0-9 .-]+)$", head.strip())
    if serovar and head.split()[:1] == ["Salmonella"]:
        add(f"Salmonella enterica subsp. enterica serovar {serovar.group(1).strip()}")
    stripped = _SET_QUALIFIERS.sub(" ", head)
    add(stripped)
    words = [w for w in re.split(r"[\s/]+", stripped) if w]
    # A binomial anywhere in the row: "Hybrid diarrheagenic Escherichia coli
    # pathotypes" is a trait of Escherichia coli.
    for first, second in zip(words, words[1:], strict=False):
        if taxon_like(first) and second[:1].islower() and second.isalpha() and len(second) > 2:
            add(f"{first} {second}")
    # Then a bare genus, the leading word first: "Animal Giardia competitors"
    # is a set of Giardia genomes, and that row exists to absorb their reads.
    for word in words:
        if taxon_like(word):
            add(word)
    # Last resort: a genus named only inside the row's parenthetical, as in
    # "Fish ciliate competitors (including fish Balantidium)".
    for word in re.split(r"[^A-Za-z]+", name):
        if taxon_like(word):
            add(word)
    return [v for v in out if v.lower() not in _NOT_TAXON_WORDS]


def _resolve_taxon(names: Sequence[str]) -> tuple[str | None, tuple[int, ...], str]:
    """The single best taxon for a row, or a verdict. See `_resolve_candidates`."""
    candidates, state = _resolve_candidates(names)
    if candidates:
        accepted, taxids = candidates[0]
        return accepted, taxids, "resolved"
    return None, (), state


def _resolve_candidates(
    names: Sequence[str],
) -> tuple[list[tuple[str, tuple[int, ...]]], str]:
    """Every taxon a row's names resolve to, best first, with a fallback verdict.

    Returning a list rather than one answer fixes two ways coverage was lost.
    A single rate-limited lookup used to abandon the row's remaining names, so
    *Cyclospora cayetanensis* was filed as having no NCBI match during a busy
    build. And a name that resolved to a taxon with nothing deposited under it
    ended the search, so "Echinococcus granulosus sensu lato" — a real taxon
    with no assembly of its own — never fell back to *Echinococcus
    granulosus*, which has one. The caller now walks the list until something
    is actually installable.
    """
    out: list[tuple[str, tuple[int, ...]]] = []
    seen_taxids: set[tuple[int, ...]] = set()
    tried: list[str] = []
    lookup_failed = False
    ambiguous: tuple[int, ...] = ()
    for requested in names:
        for name in _name_variants(requested):
            if name in tried:
                continue
            tried.append(name)
            # A network or service failure is not a statement about the
            # organism, and under a parallel build these are rate limits.
            # Backed off, retried, and on final failure noted and skipped —
            # never recorded as "this taxon does not exist".
            doc = None
            for pause in (0.0, 2.0, 6.0):
                if pause:
                    time.sleep(pause)
                try:
                    doc = _get_json(f"{_DATASETS}/taxonomy/taxon/{_quote(name)}")
                    break
                except ReferenceFetchError:
                    continue
            if doc is None:
                lookup_failed = True
                continue
            accepted, taxids, state = _nodes_to_taxon(doc.get("taxonomy_nodes") or [], name)
            if state == "ambiguous_name":
                ambiguous = ambiguous or taxids
                continue
            if state != "resolved" or not accepted:
                continue
            if accepted.strip().lower() in _FORBIDDEN_TAXA:
                # The host, or a rank so high it means nothing. Every plainer
                # spelling can only be higher, so this row is done.
                return out, "resolved_to_host_or_high_rank"
            if taxids not in seen_taxids:
                seen_taxids.add(taxids)
                out.append((accepted, taxids))
    if out:
        return out, "resolved"
    if ambiguous:
        return [], "ambiguous_name"
    return [], "lookup_unavailable" if lookup_failed else "unresolved_taxonomy"


def _nodes_to_taxon(
    nodes: Sequence[Mapping[str, Any]], queried: str = ""
) -> tuple[str | None, tuple[int, ...], str]:
    """One NCBI taxonomy response to a single accepted taxon, or a verdict.

    Two kinds of collision are decided here rather than sent to review,
    because in both the answer is determined and refusing it cost real
    coverage: *Plasmodium* and *Leishmania* went unsearched because each name
    returns a genus and a subgenus inside it, and *Giardia* because the name
    is also a retired synonym of a freshwater snail.

    * A name that is the *current* name of one candidate and a synonym of the
      others belongs to that candidate.
    * When the survivors are nested, the outermost is the answer: a row asking
      for "Plasmodium spp." wants the genus that contains the subgenus.

    Anything still tied is a genuine homonym and stays a review item.
    """
    found: dict[int, str] = {}
    lineages: dict[int, tuple[int, ...]] = {}
    for node in nodes:
        if node.get("errors"):
            continue
        tax = node.get("taxonomy") or {}
        taxid = tax.get("tax_id")
        sci = tax.get("organism_name") or (
            (tax.get("current_scientific_name") or {}).get("name")
        )
        if taxid and sci:
            found[int(taxid)] = str(sci)
            lineages[int(taxid)] = tuple(int(x) for x in (tax.get("lineage") or []) if str(x).isdigit())
    if len(found) > 1 and queried:
        exact = {t: n for t, n in found.items() if n.strip().lower() == queried.strip().lower()}
        if exact:
            found = exact
    if len(found) > 1:
        outermost = {
            t for t in found
            if not any(other != t and t in lineages.get(other, ()) for other in found)
            or all(other == t or other in lineages.get(t, ()) for other in found)
        }
        ancestors = {t for t in found if all(other == t or t in lineages.get(other, ()) for other in found)}
        if len(ancestors) == 1:
            found = {t: found[t] for t in ancestors}
        elif len(outermost) == 1:
            found = {t: found[t] for t in outermost}
    if len(found) == 1:
        taxid, sci = next(iter(found.items()))
        return sci, (taxid,), "resolved"
    if len(found) > 1:
        return None, tuple(sorted(found)), "ambiguous_name"
    return None, (), "unresolved_taxonomy"


def _best_assembly(taxid: int) -> tuple[dict[str, Any] | None, tuple[str, ...]]:
    """Pick one assembly for a taxon, deterministically.

    Preference order: RefSeq reference genome, then representative, then the
    most complete assembly, then the largest — with accession as the final
    tie-break so the choice is reproducible across builds.
    """
    url = (
        f"{_DATASETS}/genome/taxon/{taxid}/dataset_report"
        "?filters.assembly_version=current&page_size=40"
        "&filters.exclude_atypical=true"
    )
    try:
        doc = _get_json(url)
    except ReferenceFetchError:
        return None, ()
    reports = doc.get("reports") or []
    if not reports:
        return None, ()

    def key(report: Mapping[str, Any]) -> tuple[int, int, int, int, str]:
        info = report.get("assembly_info") or {}
        stats = report.get("assembly_stats") or {}
        category = str(info.get("refseq_category", "")).lower()
        level = str(info.get("assembly_level", ""))
        is_refseq = 0 if str(report.get("source_database", "")).endswith("REFSEQ") else 1
        try:
            length = int(stats.get("total_sequence_length") or 0)
        except (TypeError, ValueError):
            length = 0
        return (
            _CATEGORY_RANK.get(category, 9),
            _LEVEL_RANK.get(level, 9),
            is_refseq,
            -length,
            str(report.get("accession", "")),
        )

    ordered = sorted(reports, key=key)
    considered = tuple(str(r.get("accession", "")) for r in ordered[:8])
    return ordered[0], considered


def _pinned_assembly(accession: str) -> tuple[dict[str, Any] | None, tuple[str, ...]]:
    """Fetch one named assembly. A pin is a decision, not a preference."""
    url = f"{_DATASETS}/genome/accession/{_quote(accession)}/dataset_report"
    try:
        doc = _get_json(url)
    except ReferenceFetchError:
        return None, ()
    reports = doc.get("reports") or []
    if not reports:
        return None, ()
    return reports[0], (str(reports[0].get("accession", "")),)


def _marker_accessions(name: str, group: str, *, limit: int = 4) -> tuple[str, ...]:
    """Find curated organelle or ribosomal marker records for one organism.

    Marker evidence is a legitimate, separately-tracked route. It is *not* a
    substitute for a genome: everything installed this way is capped at
    `marker_signal` no matter how persuasive it looks.
    """
    found: list[str] = []
    for term in _MARKER_QUERIES.get(group, ("complete genome",)):
        query = f'"{name}"[Organism] AND {term}[Title] AND 200:200000[SLEN]'
        url = (
            f"{_EUTILS}/esearch.fcgi?db=nuccore&retmode=json&retmax=3"
            f"&term={_quote(query)}"
        )
        try:
            doc = _get_json(url, timeout=60.0)
        except ReferenceFetchError:
            continue
        ids = ((doc.get("esearchresult") or {}).get("idlist")) or []
        found.extend(str(i) for i in ids)
        if len(found) >= limit:
            break
        time.sleep(0.12)  # stay inside the unauthenticated E-utilities rate
    # De-duplicate, preserving discovery order for reproducibility.
    seen: set[str] = set()
    unique = [i for i in found if not (i in seen or seen.add(i))]
    return tuple(unique[:limit])


#: How specific a catalogue row's rank is. When two rows would install the
#: same assembly, the more specific one keeps it: a named species is what a
#: reader can act on, and a genus umbrella over it adds nothing a search can
#: use.
_RANK_SPECIFICITY: Final[Mapping[str, int]] = {
    "strain": 0, "serovar": 1, "serotype": 1, "pathotype": 1, "biotype": 1,
    "subspecies": 2, "species": 3, "species_complex": 4, "group": 5, "genus": 6, "family": 7,
}


#: Below this much target-specific sequence a genome-backed row is, in
#: practice, uncallable: a handful of k-mers cannot clear a breadth or
#: region-count rule. Named at build time rather than discovered as a silent
#: negative on every report.
MIN_TARGET_INFORMATIVE_BASES: Final = 10_000


def deduplicate_references(
    catalog: Catalog,
    resolutions: MutableMapping[str, Resolution],
    *,
    progress: Callable[[str], None] | None = None,
) -> int:
    """Refuse to install one assembly under two rows. Returns rows folded.

    Two rows over the same genome are not a storage inefficiency. The
    ownership pass decides which sequence is specific to one target by
    comparing every target against every other, so an identical second copy
    makes *both* rows non-specific: each loses the sequence that identified
    it, and each is then reported as unresolvable. Twelve targets were lost
    that way when genus-level rows were first installed, among them *Giardia
    duodenalis*, *Cryptosporidium parvum* and *Aspergillus flavus* — the
    genus row's best assembly was, in each case, literally the species row's
    assembly.

    The more specific row keeps the reference; the other records which row
    covers it, so the report can say it was searched and at what rank.
    """
    seeds = {s.seed_id: s for s in catalog.seeds}

    #: Which row a reader needs named when several share one genome. The
    #: recognised human pathogen, ahead of an animal or environmental
    #: relative kept only to absorb its reads: *Cyclospora cayetanensis* is
    #: the notifiable foodborne parasite, and an alphabetical tie-break gave
    #: its genome to *Cyclospora ashfordi*, a lineage nobody screens for.
    clinical = {
        "established_enteric": 0, "toxin_or_pathotype_dependent": 1, "conditional_opportunist": 2,
        "rare_enteric": 3, "extraintestinal_watch": 4, "uncertain_enteric_role": 5,
        "background_or_decoy": 6,
    }

    def key(seed_id: str) -> tuple[int, int, int, int, int, str]:
        seed = seeds[seed_id]
        accepted = (resolutions[seed_id].accepted_taxon or "").strip()
        # A row that resolved to a *named species* outranks one that fell back
        # to the genus. *Cyclospora ashfordi* and *henanensis* are lineage
        # names NCBI does not carry, so both resolve to the genus, while
        # *Cyclospora cayetanensis* — the notifiable foodborne parasite — is a
        # real taxon; without this the genome was filed under a lineage
        # nobody screens for.
        resolved_to_species = 0 if len(accepted.split()) >= 2 else 1
        return (
            _RANK_SPECIFICITY.get(seed.taxonomic_resolution, 9),
            1 if seed.expansion_state == "manual_review_required" else 0,
            resolved_to_species,
            clinical.get(seed.interpretation_class, 9),
            # A row written as a variant of another ("sensu lato", "lineage
            # B") yields to the plainly named organism.
            1 if _SET_QUALIFIERS.search(seed.requested_name or "") else 0,
            seed_id,
        )

    groups: dict[str, list[str]] = {}
    for seed_id, res in resolutions.items():
        if res.route == "genome" and res.accession:
            groups.setdefault(f"genome:{res.accession}", []).append(seed_id)
        elif res.route == "marker" and res.marker_accessions:
            groups.setdefault("marker:" + ",".join(sorted(res.marker_accessions)), []).append(seed_id)

    folded = 0
    for asset, members in sorted(groups.items()):
        if len(members) < 2:
            continue
        keeper, *rest = sorted(members, key=key)
        for seed_id in rest:
            res = resolutions[seed_id]
            res.route = "none"
            res.resolution_state = "covered_by_relative"
            res.gap_reason = "covered_by_relative"
            res.covered_by = keeper
            res.accession = None
            res.marker_accessions = ()
            res.note = (
                f"The reference for this row is the same {asset.split(':', 1)[0]} as "
                f"{seeds[keeper].display_name!r} ({asset.split(':', 1)[1]}). Installing it twice "
                "would leave both rows with no sequence of their own, so the reads are searched "
                f"against it once, under {keeper}, and this row is reported at the rank the two share."
            )
            folded += 1
            if progress is not None:
                progress(f"folded {seed_id} into {keeper} (same {asset.split(':', 1)[0]})")
    return folded


def _install_report(res: Resolution, report: Mapping[str, Any], *, max_asset_bytes: int) -> bool:
    """Record one assembly on a resolution. True when it is the genome route."""
    stats = report.get("assembly_stats") or {}
    info = report.get("assembly_info") or {}
    try:
        size = int(stats.get("total_sequence_length") or 0)
    except (TypeError, ValueError):
        size = 0
    res.accession = str(report.get("accession") or "") or None
    res.assembly_level = str(info.get("assembly_level") or "") or None
    res.assembly_bytes = size
    if size and size <= max_asset_bytes:
        res.route = "genome"
        return True
    if size:
        # Over budget: the nuclear genome exists but is not in this bundle, and
        # the row falls to organelle and ribosomal markers.
        res.note = (
            f"Nuclear assembly {res.accession} is {size / 1e6:.0f} Mb, above this bundle's "
            f"{max_asset_bytes / 1e6:.0f} Mb per-asset budget; organelle and ribosomal "
            "markers were installed instead."
        )
    return False


def _install_markers(res: Resolution, name: str, group: str, *, max_asset_bytes: int) -> bool:
    """Try the marker route for one taxon. True when markers were found."""
    markers = _marker_accessions(name, group)
    if not markers:
        return False
    res.marker_accessions = markers
    res.route = "marker"
    if res.accession and res.assembly_bytes > max_asset_bytes:
        res.gap_reason = "nuclear_genome_over_asset_budget"
    else:
        res.gap_reason = "no_whole_genome_assembly_available"
        res.accession = None
    return True


def resolve_catalog(
    catalog: Catalog,
    *,
    max_asset_bytes: int = 600_000_000,
    workers: int = 8,
    cache: Path | None = None,
    progress: Callable[[str], None] | None = None,
) -> dict[str, Resolution]:
    """Resolve every seed in the catalog. Never fails for one bad seed.

    Malformed build inputs fail; an organism that simply has no usable public
    reference produces an explicit unresolved record and the rest of the
    catalog continues (spec §12.2).
    """
    cached: dict[str, Any] = {}
    if cache is not None and cache.is_file():
        try:
            cached = json.loads(cache.read_text())
        except (OSError, json.JSONDecodeError):
            cached = {}

    #: Verdicts worth keeping between builds. A failure is not: four of the
    #: rows this report called "no NCBI match" - Enterobius vermicularis
    #: (pinworm), Cyclospora cayetanensis, Campylobacter showae, Vibrio
    #: furnissii - resolve on the first attempt today. The lookups had failed
    #: transiently and the failure was cached as though it were a fact about
    #: the organism. Only success is cached; every gap is retried.
    def one(seed: Seed) -> Resolution:
        prior = cached.get(seed.seed_id)
        if prior and str(prior.get("route") or "none") in {"genome", "marker"}:
            return Resolution(**prior)

        res = Resolution(seed_id=seed.seed_id, requested_name=seed.requested_name)

        # A genus- or group-level row flagged for manual review used to stop
        # here, unresolved and unsearched. That left 55 rows with no reference
        # at all — Blastocystis, Cryptococcus neoformans, Histoplasma,
        # Fusarium, Plasmodium among them — reported as "not searched" when
        # the reason was that nobody had expanded the row into species. The
        # row declares its own resolution (`genus`, `species_complex`), every
        # result carries `species_resolution_limited`, and the report states
        # the rank, so searching the row's own taxon cannot become a species
        # claim. The review flag stays on the record; it no longer means the
        # organism is skipped.
        awaiting_review = seed.expansion_state == "manual_review_required"

        # A trait, complex or group row resolves against its declared
        # representative taxon first; only then do the row's own name and
        # aliases get a turn.
        names = [
            n
            for n in (seed.resolution_name, seed.requested_name, *seed.aliases)
            if n
        ]
        candidates, state = _resolve_candidates(names)
        accepted, taxids = candidates[0] if candidates else (None, ())
        res.accepted_taxon = accepted
        res.ncbi_taxids = taxids
        res.resolution_state = state
        if state == "ambiguous_name":
            res.gap_reason = "ambiguous_taxonomy_requires_review"
            res.note = (
                f"{seed.requested_name!r} matched {len(taxids)} distinct taxa; "
                "a first-match selection is not permitted."
            )
            return res
        if state == "resolved_to_host_or_high_rank":
            res.resolution_state = "no_usable_reference"
            res.gap_reason = "no_species_level_taxon_for_this_row"
            res.note = (
                f"{seed.requested_name!r} has no taxon of its own: the plainest spelling of it "
                "resolves to the host or to a kingdom-level rank, which cannot be a target."
            )
            return res
        if state == "lookup_unavailable":
            res.resolution_state = "lookup_unavailable"
            res.gap_reason = "taxonomy_lookup_unavailable"
            res.note = (
                f"NCBI taxonomy did not answer for {seed.requested_name!r}. This is a "
                "service failure, not a statement about the organism; the next build retries it."
            )
            return res
        if state != "resolved" or not taxids:
            res.gap_reason = "unresolved_taxonomy"
            res.note = (
                f"No NCBI taxonomy match for {seed.requested_name!r} or any plainer "
                f"spelling of it ({', '.join(_name_variants(seed.requested_name)[1:]) or 'none'})."
            )
            return res
        if awaiting_review:
            res.note = (
                f"Group-level row searched against its declared taxon {accepted!r}. Members are "
                "not separated: every result from this row is reported at that rank."
            ).strip()

        if seed.pinned_accession:
            report, considered = _pinned_assembly(seed.pinned_accession)
            if report is None:
                res.resolution_state = "no_usable_reference"
                res.route = "none"
                res.gap_reason = "pinned_accession_unavailable"
                res.note = (
                    f"{seed.pinned_accession} was pinned for this row but NCBI returned no "
                    "current report for it. A pin is deliberate, so no substitute is chosen."
                )
                return res
            res.note = (
                f"assembly {seed.pinned_accession} pinned by the catalogue rather than chosen "
                "as the taxon's best"
            )
            res.candidates_considered = considered
            if _install_report(res, report, max_asset_bytes=max_asset_bytes):
                return res
            if _install_markers(res, accepted or seed.requested_name, seed.group,
                                max_asset_bytes=max_asset_bytes):
                return res
            res.route = "none"
            res.resolution_state = "no_usable_reference"
            res.gap_reason = "no_usable_reference"
            return res

        # Walk the candidate taxa until one is actually installable. A taxon
        # that exists but has nothing deposited under it is not an answer:
        # "Echinococcus granulosus sensu lato" is a real NCBI taxon with no
        # assembly, and stopping there left a serious parasite unsearched
        # while *Echinococcus granulosus*, one variant plainer, has a genome.
        for accepted_name, candidate_taxids in candidates:
            report, considered = _best_assembly(candidate_taxids[0])
            res.candidates_considered = considered
            if report is not None:
                res.accepted_taxon = accepted_name
                res.ncbi_taxids = candidate_taxids
                if _install_report(res, report, max_asset_bytes=max_asset_bytes):
                    return res
            if _install_markers(res, accepted_name, seed.group, max_asset_bytes=max_asset_bytes):
                res.accepted_taxon = accepted_name
                res.ncbi_taxids = candidate_taxids
                return res

        res.route = "none"
        res.resolution_state = "no_usable_reference"
        res.gap_reason = "no_usable_reference"
        tried = " or ".join(repr(name) for name, _ in candidates) or repr(seed.requested_name)
        res.note = (
            f"No current assembly and no curated marker record for {tried}. This is a "
            "reference gap, not a negative result: nothing has been deposited to search against."
        )
        return res

    results: dict[str, Resolution] = {}
    with futures.ThreadPoolExecutor(max_workers=workers) as pool:
        jobs = {pool.submit(one, seed): seed for seed in catalog.seeds}
        for done, job in enumerate(futures.as_completed(jobs), 1):
            seed = jobs[job]
            try:
                results[seed.seed_id] = job.result()
            except Exception as exc:  # noqa: BLE001 — one seed cannot stop a build
                results[seed.seed_id] = Resolution(
                    seed_id=seed.seed_id,
                    requested_name=seed.requested_name,
                    resolution_state="resolution_failed",
                    gap_reason="resolution_error",
                    note=f"{type(exc).__name__}: {exc}",
                )
            if progress is not None and done % 25 == 0:
                progress(f"resolved {done}/{len(catalog.seeds)}")

    if cache is not None:
        cache.parent.mkdir(parents=True, exist_ok=True)
        cache.write_text(
            json.dumps({k: v.to_json() for k, v in sorted(results.items())}, indent=1)
        )
    return results


# --------------------------------------------------------------------------- #
# Acquisition
# --------------------------------------------------------------------------- #


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 22), b""):
            digest.update(block)
    return digest.hexdigest()


def _download_genome(accession: str, dest: Path) -> tuple[Path, str]:
    """Fetch one assembly's genomic FASTA by accession.version.

    The Datasets download endpoint returns a zip. We extract only the FASTA
    members and never execute anything from the archive: a reference package
    is untrusted scientific input (spec §15.3).
    """
    url = (
        f"{_DATASETS}/genome/accession/{_quote(accession)}/download"
        "?include_annotation_type=GENOME_FASTA&filename=pkg.zip"
    )
    tmp = dest.with_suffix(".zip.tmp")
    net.download_to(url, tmp, timeout=1800.0)
    dest.parent.mkdir(parents=True, exist_ok=True)
    wrote = 0
    try:
        with zipfile.ZipFile(tmp) as archive, dest.open("wb") as out:
            members = [
                m
                for m in archive.namelist()
                if m.endswith((".fna", ".fa", ".fasta")) and "/data/" in m
            ]
            if not members:
                raise ReferenceFetchError(
                    f"{accession}: download contained no genomic FASTA"
                )
            for member in sorted(members):
                with archive.open(member) as handle:
                    wrote += out.write(handle.read())
    finally:
        tmp.unlink(missing_ok=True)
    if not wrote:
        dest.unlink(missing_ok=True)
        raise ReferenceFetchError(f"{accession}: empty genomic FASTA")
    return dest, _sha256(dest)


def _download_markers(ids: Sequence[str], dest: Path) -> tuple[Path, str]:
    """Fetch marker nucleotide records by UID via E-utilities efetch."""
    url = (
        f"{_EUTILS}/efetch.fcgi?db=nuccore&rettype=fasta&retmode=text"
        f"&id={_quote(','.join(ids))}"
    )
    response = net.get(url, timeout=180.0, attempts=4)
    body = response.body
    if not body.lstrip().startswith(b">"):
        raise ReferenceFetchError(f"markers {ids}: response was not FASTA")
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(body)
    return dest, _sha256(dest)


@dataclass
class Asset:
    """One downloaded, checksummed reference file."""

    seed_id: str
    path: Path
    sha256: str
    sequence_type: str
    accession_version: str
    source_url: str
    retrieved_at: str
    bases: int = 0
    license_id: str = "ncbi-public-domain-with-possible-upstream-rights"

    def to_json(self) -> dict[str, Any]:
        return {
            "seed_id": self.seed_id,
            "path": str(self.path),
            "sha256": self.sha256,
            "sequence_type": self.sequence_type,
            "accession_version": self.accession_version,
            "source_url": self.source_url,
            "retrieved_at": self.retrieved_at,
            "bases": self.bases,
            "license_id": self.license_id,
        }


def acquire(
    resolutions: Mapping[str, Resolution],
    dest: Path,
    *,
    workers: int = 6,
    progress: Callable[[str], None] | None = None,
) -> dict[str, Asset]:
    """Download every resolved reference. Failures become gaps, not crashes."""
    dest.mkdir(parents=True, exist_ok=True)
    wanted = [r for r in resolutions.values() if r.route in {"genome", "marker"}]

    def one(res: Resolution) -> Asset | None:
        target_dir = dest / res.seed_id
        stamp = target_dir / "asset.json"
        if stamp.is_file():
            try:
                prior = json.loads(stamp.read_text())
                path = Path(prior["path"])
                if path.is_file() and _sha256(path) == prior["sha256"]:
                    return Asset(**{**prior, "path": path})
                # Checksum mismatch invalidates this asset alone.
                if progress is not None:
                    progress(f"{res.seed_id}: checksum mismatch, re-downloading")
            except (OSError, KeyError, json.JSONDecodeError):
                pass
        target_dir.mkdir(parents=True, exist_ok=True)
        if res.route == "genome" and res.accession:
            path, digest = _download_genome(res.accession, target_dir / "genome.fna")
            asset = Asset(
                seed_id=res.seed_id,
                path=path,
                sha256=digest,
                sequence_type="genome",
                accession_version=res.accession,
                source_url=f"{_DATASETS}/genome/accession/{res.accession}/download",
                retrieved_at=ReferenceLock.now(),
            )
        elif res.route == "marker" and res.marker_accessions:
            path, digest = _download_markers(
                res.marker_accessions, target_dir / "markers.fna"
            )
            asset = Asset(
                seed_id=res.seed_id,
                path=path,
                sha256=digest,
                sequence_type="marker",
                accession_version=",".join(res.marker_accessions),
                source_url=f"{_EUTILS}/efetch.fcgi?db=nuccore",
                retrieved_at=ReferenceLock.now(),
            )
        else:
            return None
        asset.bases = sum(len(seq) for _, seq in iter_fasta(asset.path))
        stamp.write_text(json.dumps(asset.to_json(), indent=1))
        return asset

    assets: dict[str, Asset] = {}
    with futures.ThreadPoolExecutor(max_workers=workers) as pool:
        jobs = {pool.submit(one, r): r for r in wanted}
        for done, job in enumerate(futures.as_completed(jobs), 1):
            res = jobs[job]
            try:
                asset = job.result()
            except Exception as exc:  # noqa: BLE001
                res.gap_reason = "download_failed"
                res.note = f"{type(exc).__name__}: {exc}"
                res.route = "none"
                asset = None
            if asset is not None:
                assets[res.seed_id] = asset
            if progress is not None and done % 20 == 0:
                progress(f"acquired {done}/{len(wanted)}")
    return assets


# --------------------------------------------------------------------------- #
# Bundle build
# --------------------------------------------------------------------------- #


@dataclass
class BundleManifest:
    """The frozen installed manifest a sample run consumes."""

    bundle_id: str
    catalog_version: str
    entries: dict[str, dict[str, Any]] = field(default_factory=dict)
    lock: dict[str, Any] = field(default_factory=dict)
    decoys: tuple[str, ...] = ()

    def to_json(self) -> dict[str, Any]:
        return {
            "bundle_id": self.bundle_id,
            "catalog_version": self.catalog_version,
            "entries": self.entries,
            "lock": self.lock,
            "decoys": list(self.decoys),
        }

    @classmethod
    def load(cls, path: Path) -> BundleManifest:
        doc = json.loads(Path(path).read_text())
        return cls(
            bundle_id=str(doc["bundle_id"]),
            catalog_version=str(doc["catalog_version"]),
            entries=dict(doc.get("entries") or {}),
            lock=dict(doc.get("lock") or {}),
            decoys=tuple(doc.get("decoys") or ()),
        )


def _bowtie2_version(binary: str = "bowtie2") -> str:
    exe = shutil.which(binary)
    if exe is None:
        return "not installed"
    try:
        out = subprocess.run(  # noqa: S603 - fixed argv, no shell
            [exe, "--version"], capture_output=True, text=True, timeout=60, check=False
        )
        return out.stdout.splitlines()[0].strip() if out.stdout else "unknown"
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def build_bundle(
    catalog: Catalog,
    resolutions: Mapping[str, Resolution],
    assets: Mapping[str, Asset],
    *,
    out_dir: Path,
    decoy_fastas: Sequence[Path] = (),
    threads: int = 8,  # noqa: ARG001 - recorded for the competitive index step
    bowtie2_build: str = "bowtie2-build",
    progress: Callable[[str], None] | None = None,
    build_aligner_index: bool = True,
    decoy_cache: Path | None = None,
) -> BundleManifest:
    """Compute informative regions competitively and write the bundle.

    Ordering is the whole point. Every asset's k-mers are collected *first*,
    then ownership is decided across the pooled set, and only then is any
    target credited with discriminating sequence. Deciding target by target
    would credit whichever genome was processed first with sequence it shares
    with a relative — the exact failure mode behind contaminated-reference
    worm calls (spec §1, ParaRef).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    say = progress or (lambda _m: None)

    ordered = sorted(assets)
    if not ordered:
        raise PathogenError("no reference assets were acquired; nothing to build")
    if len(ordered) >= SHARED:
        raise PathogenError(
            f"{len(ordered)} targets exceeds the {SHARED - 1}-target index width"
        )
    index_of = {seed_id: i for i, seed_id in enumerate(ordered)}

    # ---- pass 1: collect (kmer, owner, bin) triples ---------------------- #
    kmer_blocks: list[np.ndarray] = []
    owner_blocks: list[np.ndarray] = []
    bin_blocks: list[np.ndarray] = []
    per_target_bases: dict[str, int] = {}

    def collect(path: Path, owner: int) -> None:
        offset = 0
        for _header, seq in iter_fasta(path):
            if len(seq) < 200:
                offset += len(seq)
                continue
            kmers, positions = select_kmers(seq, positions=True)
            if kmers.size:
                kmer_blocks.append(kmers)
                owner_blocks.append(np.full(kmers.size, owner, dtype=np.uint16))
                assert positions is not None
                bin_blocks.append(
                    ((positions + offset) // BIN_BASES).astype(np.int64)
                )
            offset += len(seq)

    for n, seed_id in enumerate(ordered, 1):
        asset = assets[seed_id]
        collect(asset.path, index_of[seed_id])
        per_target_bases[seed_id] = asset.bases
        if n % 50 == 0:
            say(f"hashed {n}/{len(ordered)} references")

    kmers = np.concatenate(kmer_blocks)
    owners = np.concatenate(owner_blocks)
    bins = np.concatenate(bin_blocks)
    # Emptied rather than deleted: `collect` closes over these names, and
    # `concatenate` has already copied, so clearing releases the blocks
    # without leaving the closure referring to unbound names.
    kmer_blocks.clear()
    owner_blocks.clear()
    bin_blocks.clear()
    say(f"{kmers.size:,} sampled reference k-mers; resolving ownership")

    # Decoys — the host, food animals, laboratory vectors and near-neighbour
    # background — exist only to *remove* sequence from every target's
    # informative set, so they need no identity and no bins of their own.
    # Holding them in a separate sorted set rather than concatenating them
    # into the main table roughly halves peak memory on a multi-gigabase
    # bundle, which matters because these genomes are the largest ones here.
    decoy_kmers = _decoy_kmer_set(decoy_fastas, say, cache=decoy_cache)

    # ---- pass 2: decide ownership across the pooled set ------------------- #
    # Introsort, not a stable sort: on ~10^9 k-mers NumPy's stable argsort
    # (timsort) is several times slower, and nothing below depends on the
    # order *within* a run of identical k-mers — ownership is an all-equal
    # test over the run and the bin is taken as the run minimum, both of
    # which are order-free. Output is therefore identical and deterministic.
    order = np.argsort(kmers, kind="quicksort")
    kmers = kmers[order]
    owners = owners[order]
    bins = bins[order]
    # The permutation is 8 bytes per k-mer — as large as the k-mer array
    # itself — and nothing below needs it. Holding it through the grouping
    # pass would add several gigabytes to peak usage for no reason.
    del order

    # Group boundaries of identical k-mers. Everything below is vectorised:
    # at ~10^9 sampled k-mers a Python loop over groups would take days.
    boundary = np.concatenate(([True], kmers[1:] != kmers[:-1]))
    starts = np.flatnonzero(boundary)
    group_of = np.cumsum(boundary) - 1
    del boundary
    first_owner = owners[starts]
    # A k-mer is target-specific only when every one of its occurrences has
    # the same owner. `logical_and.reduceat` gives that per group in one pass.
    matches_first = owners == first_owner[group_of]
    same = np.logical_and.reduceat(matches_first, starts)
    resolved_owner = np.where(same, first_owner, np.uint16(SHARED)).astype(np.uint16)
    del matches_first

    unique_kmers = kmers[starts]
    # A target-specific k-mer's occurrences all lie in one genome, and
    # `collect` emits them in ascending position, so the run minimum is the
    # first occurrence — the same bin a stable sort would have surfaced.
    # (Bins of shared k-mers are never read.)
    unique_bins = np.minimum.reduceat(bins, starts)

    # Anything the host or a food, vector or background genome also carries is
    # not target-specific, no matter how clean its within-genus uniqueness
    # looks. This is the step that stops a bacterial contig inside a worm
    # assembly from becoming worm evidence.
    if decoy_kmers.size:
        pos = np.searchsorted(decoy_kmers, unique_kmers)
        pos = np.minimum(pos, decoy_kmers.size - 1)
        in_decoy = decoy_kmers[pos] == unique_kmers
        masked = int((in_decoy & (resolved_owner != SHARED)).sum())
        resolved_owner = np.where(in_decoy, np.uint16(SHARED), resolved_owner).astype(
            np.uint16
        )
        say(f"decoys masked {masked:,} otherwise target-specific k-mers")
        del pos, in_decoy

    say(
        f"{unique_kmers.size:,} distinct k-mers, "
        f"{int((resolved_owner != SHARED).sum()):,} target-specific"
    )

    # ---- informative bins, per target ------------------------------------- #
    # Sort the target-specific k-mers by (owner, bin) once, then slice each
    # target's block with searchsorted rather than rescanning per target.
    informative_bins: dict[int, np.ndarray] = {}
    informative_kmers: dict[str, int] = dict.fromkeys(ordered, 0)
    specific = resolved_owner != SHARED
    spec_owner = resolved_owner[specific].astype(np.int64)
    spec_bin = unique_bins[specific]
    order2 = np.lexsort((spec_bin, spec_owner))
    spec_owner = spec_owner[order2]
    spec_bin = spec_bin[order2]
    lo = np.searchsorted(spec_owner, np.arange(len(ordered)), side="left")
    hi = np.searchsorted(spec_owner, np.arange(len(ordered)), side="right")
    for seed_id, owner in index_of.items():
        mine = spec_bin[lo[owner] : hi[owner]]
        informative_kmers[seed_id] = int(mine.size)
        if mine.size == 0:
            informative_bins[owner] = np.empty(0, dtype=np.int64)
            continue
        values, counts = np.unique(mine, return_counts=True)
        informative_bins[owner] = values[counts >= MIN_INFORMATIVE_KMERS_PER_BIN]
    del spec_owner, spec_bin, order2

    # ---- near neighbours from measured sharing ---------------------------- #
    neighbours = _near_neighbours(owners, group_of, same, index_of, ordered)
    del group_of

    index = KmerIndex(
        kmers=unique_kmers,
        owners=resolved_owner,
        target_ids=tuple(ordered),
        informative_bins=informative_bins,
    )
    index_path = out_dir / "kmer_index.npz"
    index.save(index_path)
    say(f"wrote {index_path.name}")

    # ---- per-target FASTA for on-demand competitive indexes --------------- #
    # The competitive stage aligns only candidate reads, so it needs an index
    # over the candidates and their measured relatives — not over every
    # installed genome. Writing one normalised FASTA per target lets that
    # index be assembled on demand and cached, which turns a multi-hour
    # whole-bundle index build into a few minutes for the targets that a
    # given sample actually nominated.
    seq_dir = out_dir / "sequences"
    seq_dir.mkdir(parents=True, exist_ok=True)
    total_written = 0
    for seed_id in ordered:
        dest = seq_dir / f"{seed_id}.fna"
        if dest.is_file() and dest.stat().st_size:
            total_written += dest.stat().st_size
            continue
        with dest.open("wb") as out:
            for header, seq in iter_fasta(assets[seed_id].path):
                contig = header.split()[0] if header.split() else "contig"
                out.write(f">{seed_id}|{contig}\n".encode())
                for i in range(0, len(seq), 80):
                    out.write(seq[i : i + 80] + b"\n")
        total_written += dest.stat().st_size
    say(f"wrote {len(ordered)} per-target FASTAs ({total_written / 1e9:.2f} GB)")

    aligner = _bowtie2_version() if shutil.which(bowtie2_build) else (
        "bowtie2-build not installed"
    )
    if build_aligner_index and shutil.which(bowtie2_build) is None:
        say("bowtie2-build missing: competitive confirmation will be unavailable")

    # ---- guard: nothing installed may end up unidentifiable --------------- #
    # A target with its own installed genome that ends with no k-mer of its own
    # has been masked by something else in this bundle, and it will be reported
    # as unresolvable for every sample. That happened silently once, to
    # *Giardia duodenalis* and eleven others, when genus-level rows duplicated
    # species assemblies. It is a build error, not a result, so the build says
    # so rather than writing a manifest that quietly drops an organism.
    blinded = sorted(
        seed_id for seed_id in ordered
        if assets[seed_id].sequence_type == "genome"
        and not informative_kmers.get(seed_id, 0)
    )
    if blinded:
        say(
            f"WARNING: {len(blinded)} installed target(s) have no sequence of their own after "
            "masking and cannot be identified from this bundle: " + ", ".join(blinded)
        )
    # Not zero, but too little to call. A genome-backed target holding a few
    # dozen specific k-mers out of millions of reference bases is one
    # near-identical neighbour away from silence, and the build should name it
    # while a human is still watching.
    thin = sorted(
        (informative_kmers.get(seed_id, 0), seed_id)
        for seed_id in ordered
        if assets[seed_id].sequence_type == "genome"
        and 0 < informative_kmers.get(seed_id, 0) * SCALED < MIN_TARGET_INFORMATIVE_BASES
    )
    if thin:
        say(
            f"WARNING: {len(thin)} installed target(s) hold under "
            f"{MIN_TARGET_INFORMATIVE_BASES:,} bases of their own sequence and will rarely be "
            "callable: " + ", ".join(f"{s} ({k * SCALED:,} bp)" for k, s in thin[:12])
        )

    # ---- manifest --------------------------------------------------------- #
    total_bases = sum(per_target_bases.values())
    fingerprint = hashlib.sha256(
        json.dumps(
            {
                "catalog": catalog.version,
                "k": index.k,
                "scaled": SCALED,
                "bin": BIN_BASES,
                "assets": {s: assets[s].sha256 for s in ordered},
            },
            sort_keys=True,
        ).encode()
    ).hexdigest()[:16]
    bundle_id = f"pathogens-dna-{fingerprint}"

    entries: dict[str, dict[str, Any]] = {}
    coverage_counts: dict[str, int] = {}
    for seed in catalog.seeds:
        res = resolutions.get(seed.seed_id)
        asset = assets.get(seed.seed_id)
        if asset is not None and res is not None:
            if asset.sequence_type == "genome":
                status = "genome_supported"
            else:
                status = (
                    "organelle_only"
                    if "mitochondrion" in asset.accession_version.lower()
                    else "marker_only"
                )
            # Sequence that shares every k-mer with a relative is installed but
            # cannot resolve the target. Say so rather than implying it can.
            if informative_kmers.get(seed.seed_id, 0) == 0:
                status = "unresolved_taxonomy"
            entry: dict[str, Any] = {
                "reference_status": status,
                "ncbi_taxids": list(res.ncbi_taxids),
                "accepted_taxon": res.accepted_taxon,
                "reference_accessions": (
                    [asset.accession_version] if asset.sequence_type == "genome" else []
                ),
                "marker_accessions": (
                    asset.accession_version.split(",")
                    if asset.sequence_type == "marker"
                    else []
                ),
                "near_neighbor_target_ids": list(neighbours.get(seed.seed_id, ())),
                "reference_bases": asset.bases,
                "informative_bases": informative_kmers.get(seed.seed_id, 0) * SCALED,
                "informative_kmers": informative_kmers.get(seed.seed_id, 0),
                "index_owner": index_of[seed.seed_id],
                "asset_sha256": asset.sha256,
                "sequence_call_profile_id": "research_dna_v1",
            }
            if status == "unresolved_taxonomy":
                entry["reference_gap_reason"] = (
                    "installed sequence carries no k-mer absent from its relatives, "
                    "the host and the decoys, so this target cannot be resolved "
                    "from it"
                )
            elif res.gap_reason:
                entry["reference_gap_reason"] = res.gap_reason
        else:
            if res is not None and res.resolution_state == "covered_by_relative":
                status = "covered_by_relative"
            elif res is None or res.resolution_state in {
                "unresolved_taxonomy",
                "ambiguous_name",
                "manual_review_required",
                "lookup_unavailable",
            }:
                status = "unresolved_taxonomy"
            else:
                status = "no_usable_reference"
            entry = {
                "reference_status": status,
                "ncbi_taxids": list(res.ncbi_taxids) if res else [],
                "accepted_taxon": res.accepted_taxon if res else None,
                "reference_accessions": [],
                "marker_accessions": [],
                "near_neighbor_target_ids": [],
                "reference_bases": 0,
                "informative_bases": 0,
                "reference_gap_reason": (
                    (res.gap_reason if res else None) or "not_in_this_bundle"
                ),
                "covered_by_target_id": (res.covered_by if res else None),
                "resolution_note": res.note if res else "",
                "sequence_call_profile_id": "research_dna_v1",
            }
        entries[seed.seed_id] = entry
        coverage_counts[status] = coverage_counts.get(status, 0) + 1

    lock = ReferenceLock(
        bundle_id=bundle_id,
        created_at=ReferenceLock.now(),
        catalog_version=catalog.version,
        taxonomy_snapshot="ncbi-datasets-v2alpha-live-resolution",
        taxonomy_sha256=hashlib.sha256(
            json.dumps(
                {s: list(r.ncbi_taxids) for s, r in sorted(resolutions.items())},
                sort_keys=True,
            ).encode()
        ).hexdigest(),
        software=(
            {"name": "bowtie2", "version": aligner, "source_commit": "", "container_digest": ""},
            {"name": "openbiota.pathogens.kmers", "version": f"k{index.k}-scaled{SCALED}",
             "source_commit": "", "container_digest": ""},
        ),
        index_parameters={
            "k": index.k,
            "scaled": SCALED,
            "bin_bases": BIN_BASES,
            "min_informative_kmers_per_bin": MIN_INFORMATIVE_KMERS_PER_BIN,
            "decoys": [Path(d).name for d in decoy_fastas],
        },
        validation_profile_ids=("research_dna_v1",),
        build_command_log_sha256=fingerprint,
        asset_count=len(ordered),
        total_reference_bases=total_bases,
        target_coverage=coverage_counts,
    )

    manifest = BundleManifest(
        bundle_id=bundle_id,
        catalog_version=catalog.version,
        entries=entries,
        lock=lock.to_json(),
        decoys=tuple(Path(d).name for d in decoy_fastas),
    )
    (out_dir / "manifest.json").write_text(json.dumps(manifest.to_json(), indent=1))
    (out_dir / "resolutions.json").write_text(
        json.dumps({k: v.to_json() for k, v in sorted(resolutions.items())}, indent=1)
    )
    say(f"bundle {bundle_id}: " + ", ".join(
        f"{k}={v}" for k, v in sorted(coverage_counts.items())
    ))
    return manifest


def _decoy_cache_key(paths: Sequence[Path]) -> str:
    """Identity of a decoy set: which files, how big, when last written."""
    parts = [f"{p.resolve()}|{p.stat().st_size}|{int(p.stat().st_mtime)}" for p in paths]
    return hashlib.sha256("\n".join(parts).encode()).hexdigest()[:16]


def _decoy_kmer_set(
    decoy_fastas: Sequence[Path],
    say: Callable[[str], None],
    *,
    cache: Path | None = None,
) -> np.ndarray:
    """Sorted k-mer set for the host and background genomes.

    Sorted, not de-duplicated: the only consumer is a `searchsorted`
    membership test, which is indifferent to repeats, and an explicit
    de-duplication pass is what we must avoid — on a ~10 Gbp decoy set
    `np.unique` reaches for a hash table (NumPy ≥ 2.3) that takes tens of
    gigabytes and the better part of an hour over a billion values. An
    in-place sort is a couple of minutes and no extra memory.

    Hashing the decoys is the slowest step of a bundle build and their
    content changes only when a genome file does, so the sorted set is
    cached beside the bundle, keyed on file identity, when `cache` is given.
    """
    present = [Path(d) for d in decoy_fastas]
    for path in list(present):
        if not path.is_file():
            say(f"decoy missing, skipped: {path}")
            present.remove(path)
    if not present:
        return np.empty(0, dtype=np.uint64)

    key = _decoy_cache_key(present)
    if cache is not None and cache.is_file():
        try:
            with np.load(cache) as z:
                if str(z["key"]) == key:
                    kmers = z["kmers"]
                    say(f"reusing {kmers.size:,} cached decoy k-mers")
                    return kmers
                say("decoy cache is for a different decoy set; rebuilding")
        except (OSError, KeyError, ValueError) as exc:
            say(f"decoy cache unreadable ({exc}); rebuilding")

    blocks: list[np.ndarray] = []
    for path in present:
        say(f"hashing decoy {path.name}")
        for _header, seq in iter_fasta(path):
            if len(seq) < 200:
                continue
            kmers, _ = select_kmers(seq)
            if kmers.size:
                blocks.append(kmers)
    if not blocks:
        return np.empty(0, dtype=np.uint64)
    joined = np.concatenate(blocks)
    del blocks
    joined.sort()
    say(f"{joined.size:,} decoy k-mers sorted")
    if cache is not None:
        cache.parent.mkdir(parents=True, exist_ok=True)
        tmp = cache.with_suffix(".tmp.npz")
        np.savez(tmp, key=np.array(key), kmers=joined)
        tmp.replace(cache)
        say(f"cached decoy k-mers to {cache.name}")
    return joined


#: Statuses whose installed sequence can compete for a read.
_USABLE_REFERENCE: Final = frozenset(
    {"genome_supported", "marker_only", "organelle_only"}
)


def competitive_members(
    manifest: BundleManifest,
    target_ids: Sequence[str] | None = None,
) -> tuple[str, ...]:
    """Which references must be present for a fair competitive alignment.

    With `target_ids`, the candidates plus the near neighbours the reference
    build measured for them. With `None`, every usable reference in the
    bundle.

    The whole-bundle set is the one worth building. Measured on a real
    stool library, a sample's candidates plus neighbours already came to
    197 of 377 references and 5.4 of 7.4 Gbp — so a per-sample index costs
    most of a whole-bundle build, and because the member set shifts from
    sample to sample the cache never hits and every sample pays it again.
    One shared index is built once and reused by every run, and it is also
    strictly the more competitive of the two.
    """
    entries = manifest.entries
    if target_ids is None:
        return tuple(
            sorted(
                t for t, e in entries.items()
                if (e or {}).get("reference_status") in _USABLE_REFERENCE
            )
        )
    wanted: set[str] = set()
    for target_id in target_ids:
        entry = entries.get(target_id) or {}
        if entry.get("reference_status") in _USABLE_REFERENCE:
            wanted.add(target_id)
        for neighbour in entry.get("near_neighbor_target_ids") or ():
            nb = entries.get(neighbour) or {}
            if nb.get("reference_status") in _USABLE_REFERENCE:
                wanted.add(str(neighbour))
    return tuple(sorted(wanted))


def build_competitive_index(
    bundle_dir: Path,
    target_ids: Sequence[str] | None,
    *,
    manifest: BundleManifest,
    include_host: Path | None = None,
    threads: int = 8,
    bowtie2_build: str = "bowtie2-build",
    progress: Callable[[str], None] | None = None,
) -> Path | None:
    """Assemble and cache a bowtie2 index the candidates must compete in.

    No read is placed on a target without its real competitors present.
    Global competition is already baked in by k-mer ownership; this supplies
    the base-level competition that ownership cannot. Pass `None` for
    `target_ids` to build over the whole bundle, which is what callers
    should normally do — see :func:`competitive_members`.

    Building is slow and the result is shared, so the work is done under a
    lock and published atomically: concurrent samples wait for one build
    rather than racing into the same directory, and a build killed halfway
    leaves no half-index behind for the next run to trust.

    Returns the index prefix, or `None` when bowtie2 is unavailable.
    """
    say = progress or (lambda _m: None)
    exe = shutil.which(bowtie2_build)
    if exe is None:
        return None

    members_dir = Path(bundle_dir) / "sequences"
    wanted = competitive_members(manifest, target_ids)
    members = sorted(p for p in (members_dir / f"{t}.fna" for t in wanted) if p.is_file())
    if not members:
        return None

    fingerprint = hashlib.sha256(
        ("|".join(p.name for p in members) + f"|host={bool(include_host)}").encode()
    ).hexdigest()[:16]
    root = Path(bundle_dir) / "competitive"
    cache_dir = root / fingerprint
    prefix = cache_dir / "index"

    def ready() -> bool:
        return (cache_dir / "index.1.bt2").is_file() or (
            cache_dir / "index.1.bt2l"
        ).is_file()

    if ready():
        say(f"reusing competitive index for {len(members)} references")
        return prefix

    root.mkdir(parents=True, exist_ok=True)
    lock_path = root / f"{fingerprint}.lock"
    with lock_path.open("w") as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError:
            say("another run is building the competitive index; waiting for it")
            fcntl.flock(lock, fcntl.LOCK_EX)
        # Whoever held the lock may have finished the job for us.
        if ready():
            say(f"reusing competitive index for {len(members)} references")
            return prefix

        staging = root / f"{fingerprint}.building"
        shutil.rmtree(staging, ignore_errors=True)
        staging.mkdir(parents=True)
        try:
            combined = staging / "refs.fna"
            with combined.open("wb") as out:
                for member in members:
                    # Streamed, not slurped: the largest member genomes are
                    # hundreds of megabytes and the combined file is tens of
                    # gigabytes, so there is no reason to hold either whole.
                    with member.open("rb") as src:
                        shutil.copyfileobj(src, out, 1 << 22)
                if include_host is not None and Path(include_host).is_file():
                    # The host is the single largest source of spurious
                    # pathogen support, so it competes directly.
                    for header, seq in iter_fasta(Path(include_host)):
                        contig = header.split()[0] if header.split() else "contig"
                        out.write(f">decoy.host|{contig}\n".encode())
                        for i in range(0, len(seq), 80):
                            out.write(seq[i : i + 80] + b"\n")
            say(
                f"building competitive index over {len(members)} references "
                f"({combined.stat().st_size / 1e9:.2f} GB) — "
                "this happens once and is then reused by every sample"
            )
            proc = subprocess.run(  # noqa: S603 - fixed argv, no shell
                [exe, "--threads", str(threads), "-f", str(combined),
                 str(staging / "index")],
                capture_output=True,
                text=True,
                check=False,
            )
            combined.unlink(missing_ok=True)
            if proc.returncode != 0:
                tail = (proc.stderr or "").strip().splitlines()[-1:]
                say(f"bowtie2-build failed (rc={proc.returncode}): {' '.join(tail)}")
                return None
            shutil.rmtree(cache_dir, ignore_errors=True)
            staging.rename(cache_dir)
        finally:
            shutil.rmtree(staging, ignore_errors=True)
    say("competitive index ready")
    return prefix


def _near_neighbours(
    owners: np.ndarray,
    group_of: np.ndarray,
    same: np.ndarray,
    index_of: Mapping[str, int],
    ordered: Sequence[str],
    *,
    top: int = 6,
    sample_groups: int = 400_000,
) -> dict[str, tuple[str, ...]]:
    """Derive each target's closest relatives from measured k-mer sharing.

    These are not curated guesses: two references are neighbours here because
    they demonstrably share sequence in *this* bundle, which is exactly the
    relationship that makes a species call hard. Within-genus similarity
    alone would miss a bacterial contig sitting inside a worm assembly; this
    finds it, because the sharing is measured rather than assumed.

    A bundle can contain a billion shared k-mers, so the co-occurrence matrix
    is estimated from a deterministic random sample of shared groups. The
    ranking it produces is stable; the absolute counts are not used.
    """
    n = len(ordered)
    shared_groups = np.flatnonzero(~same)
    if shared_groups.size == 0:
        return dict.fromkeys(index_of, ())
    rng = np.random.default_rng(20260907)
    if shared_groups.size > sample_groups:
        shared_groups = np.sort(
            rng.choice(shared_groups, size=sample_groups, replace=False)
        )

    # Elements belonging to the sampled shared groups.
    keep = np.isin(group_of, shared_groups)
    g = group_of[keep]
    o = owners[keep].astype(np.int64)
    real = o != SHARED
    g, o = g[real], o[real]
    if g.size == 0:
        return dict.fromkeys(index_of, ())

    # Unique (group, owner) pairs, then all ordered owner pairs per group.
    pairs = np.unique(np.stack((g, o), axis=1), axis=0)
    g, o = pairs[:, 0], pairs[:, 1]
    bounds = np.flatnonzero(np.concatenate(([True], g[1:] != g[:-1])))
    ends = np.concatenate((bounds[1:], [g.size]))
    counts = np.zeros((n, n), dtype=np.int32)
    for start, stop in zip(bounds.tolist(), ends.tolist(), strict=True):
        block = o[start:stop]
        if block.size < 2 or block.size > 32:
            continue  # a k-mer in 30+ references says nothing about kinship
        counts[np.repeat(block, block.size), np.tile(block, block.size)] += 1
    np.fill_diagonal(counts, 0)

    out: dict[str, tuple[str, ...]] = {}
    for seed_id, owner in index_of.items():
        row = counts[owner]
        if not row.any():
            out[seed_id] = ()
            continue
        best = np.argsort(-row, kind="stable")[:top]
        out[seed_id] = tuple(ordered[int(b)] for b in best if row[int(b)] > 0)
    return out


# --------------------------------------------------------------------------- #
# Audit
# --------------------------------------------------------------------------- #


def audit_bundle(bundle_dir: Path) -> dict[str, Any]:
    """Re-derive coverage from artifacts and fail any unsupported claim.

    The manifest is a claim; the index and the asset checksums are the
    evidence. This recomputes the claim and reports every disagreement, which
    is what makes "covered" mean something (spec §12.2 `refs audit`).
    """
    bundle_dir = Path(bundle_dir)
    manifest = BundleManifest.load(bundle_dir / "manifest.json")
    index = KmerIndex.load(bundle_dir / "kmer_index.npz")
    problems: list[str] = []

    installed = set(index.target_ids)
    for seed_id, entry in manifest.entries.items():
        status = str(entry.get("reference_status"))
        claims_sequence = status in {"genome_supported", "marker_only", "organelle_only"}
        if claims_sequence and seed_id not in installed:
            problems.append(
                f"{seed_id}: claims {status} but has no sequence in the k-mer index"
            )
        if claims_sequence and not entry.get("informative_bases"):
            problems.append(
                f"{seed_id}: claims {status} with zero target-specific bases"
            )
        if claims_sequence and not (
            entry.get("reference_accessions") or entry.get("marker_accessions")
        ):
            problems.append(f"{seed_id}: claims {status} with no accession recorded")
        if not claims_sequence and not entry.get("reference_gap_reason"):
            problems.append(f"{seed_id}: is {status} but records no gap reason")

    # Asset checksums: a mismatch invalidates that asset and everything
    # derived from it, and nothing else.
    checked = 0
    for seed_id, entry in manifest.entries.items():
        digest = entry.get("asset_sha256")
        if not digest:
            continue
        stamp = bundle_dir.parent / "assets" / seed_id / "asset.json"
        if not stamp.is_file():
            problems.append(f"{seed_id}: asset record missing from disk")
            continue
        try:
            prior = json.loads(stamp.read_text())
            path = Path(prior["path"])
        except (OSError, KeyError, json.JSONDecodeError):
            problems.append(f"{seed_id}: asset record unreadable")
            continue
        if not path.is_file():
            problems.append(f"{seed_id}: asset file missing")
        elif _sha256(path) != digest:
            problems.append(f"{seed_id}: asset checksum mismatch — bundle invalid")
        checked += 1

    coverage: dict[str, int] = {}
    for entry in manifest.entries.values():
        status = str(entry.get("reference_status"))
        coverage[status] = coverage.get(status, 0) + 1

    return {
        "bundle_id": manifest.bundle_id,
        "catalog_version": manifest.catalog_version,
        "targets": len(manifest.entries),
        "assets_checksum_verified": checked,
        "indexed_targets": len(installed),
        "coverage": coverage,
        "problems": problems,
        "ok": not problems,
    }
