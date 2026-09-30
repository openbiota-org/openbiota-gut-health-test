"""Beta-glucuronidase structural classes — spec 0.8.3 A08, source F39.

A total beta-glucuronidase count answers a question nobody asked. The
enzymes in that total are not interchangeable: Pollet et al. showed that gut
GUS enzymes fall into structural classes defined by two active-site loops,
and that the loop is what decides which glucuronides an enzyme is fast on.
The enzymes that efficiently reactivate the irinotecan metabolite SN-38 are
Loop 1 enzymes; the inhibitors developed against them do not touch the
others. So "you have a lot of GUS" is not a statement about any substrate.

The class is a property of the reference protein, not of the read. That
distinction is load-bearing: *Faecalibacterium prausnitzii* encodes both a
Loop 1 GUS and the FMN-binding No Loop enzyme whose structure was solved
from strain L2-6, so classifying by organism would misassign the single
largest contributor in a typical sample. This module therefore classifies
every reference sequence, and buckets fragments by the reference each one
actually hit.

Classification follows the published definition literally. Every reference
is aligned to *E. coli* GUS (P05804) and the number of residues occupying
each of the two loop positions is counted:

  * L1 at the 356-380 region: >15 residues is L1, 10-15 is mini-Loop 1
  * L2 at the 416-419 region: >=12 residues is L2, 9-12 is mini-Loop 2

A reference whose alignment does not cover both positions is "no coverage",
which is the paper's own category for exactly this case and is reported
rather than guessed at.

What this does **not** do is estimate a drug or hormone exposure. The class
constrains which substrates are plausible; it does not measure activity,
which depends on expression, on the enzyme's kinetics against the specific
glucuronide, and on whether that glucuronide reaches the colon at all.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

METHOD_GUS_CLASSES: Final = "ext083.gus_structural_classes/1.0"

#: E. coli beta-glucuronidase, the enzyme every position in the literature
#: is numbered against.
ANCHOR_ACCESSION: Final = "P05804"

#: The two loop positions, in anchor numbering, inclusive.
L1_WINDOW: Final[tuple[int, int]] = (356, 380)
L2_WINDOW: Final[tuple[int, int]] = (416, 419)

#: Thresholds, taken from the published definitions rather than tuned.
L1_FULL_MIN: Final = 16      # ">15 residues"
L1_MINI_MIN: Final = 10      # "10-15 residues"
L2_FULL_MIN: Final = 12      # ">=12 residues"
L2_MINI_MIN: Final = 9       # "9-12 residues"

GAP: Final = "-"

#: The entry whose references these are.
GUS_ENTRY_KEY: Final = "bglucuronidase:GUS"
GUS_ENTRY_ID: Final = "GUS"


class GusClassError(ValueError):
    """A classification that could not be trusted."""


@dataclass(frozen=True)
class GusClass:
    """One structural class, and what the structural work established."""

    code: str
    label: str
    #: What the loop does, stated as the structural work states it.
    structure: str
    #: Substrate behaviour that has actually been measured for this class.
    substrates: str
    #: The claim a reader is most likely to make from this row, refused.
    does_not_mean: str

    def to_json(self) -> dict[str, Any]:
        return {
            "code": self.code, "label": self.label, "structure": self.structure,
            "substrates": self.substrates, "does_not_mean": self.does_not_mean,
        }


CLASSES: Final[tuple[GusClass, ...]] = (
    GusClass(
        code="L1", label="Loop 1",
        structure=(
            "A loop of more than 15 residues at the active site, in the position "
            "occupied by residues 356-380 of the E. coli enzyme. E. coli GUS itself is "
            "the type example, and the Clostridium perfringens and Streptococcus "
            "agalactiae enzymes have 22- and 21-residue versions."
        ),
        substrates=(
            "These are the enzymes measured to be efficient on small-molecule "
            "glucuronides, including SN-38 glucuronide, the conjugate of the irinotecan "
            "metabolite. The selective bacterial GUS inhibitors developed to reduce "
            "irinotecan gut toxicity were designed against this class and do not "
            "inhibit the others equally."
        ),
        does_not_mean=(
            "This does not estimate any drug exposure. It is a count of genes in stool, "
            "not a measurement of enzyme activity, and nobody should change a medication "
            "on the strength of it."
        ),
    ),
    GusClass(
        code="mL1", label="mini-Loop 1",
        structure=(
            "A shorter loop, 10 to 15 residues, in the same active-site position. The "
            "Bacteroides fragilis enzyme, whose structure is solved as PDB 3CMG, has a "
            "12-residue version and is the type example. This is one of the most "
            "abundant classes in the human gut."
        ),
        substrates=(
            "The shortened loop makes these enzymes generally less efficient on the "
            "small glucuronides that Loop 1 enzymes process quickly, though the "
            "relationship is not absolute and varies between individual enzymes."
        ),
        does_not_mean=(
            "A shorter loop is not an inactive enzyme. These are abundant, functional "
            "glucuronidases; the loop changes which substrates they are fast on."
        ),
    ),
    GusClass(
        code="L2", label="Loop 2",
        structure=(
            "No Loop 1, but an insertion of 12 or more residues at a second position, "
            "the 416-419 region of the E. coli enzyme. A Bacteroides uniformis enzyme is "
            "the type example, with a 19-residue loop."
        ),
        substrates=(
            "Loop 2 enzymes frequently carry additional carbohydrate-binding domains at "
            "the C-terminus that Loop 1 enzymes lack, which is consistent with acting on "
            "larger glucuronide-containing carbohydrates rather than on small drug "
            "conjugates."
        ),
        does_not_mean=(
            "The additional domains are a structural observation. They do not establish "
            "what any of these enzymes does in your gut."
        ),
    ),
    GusClass(
        code="mL2", label="mini-Loop 2",
        structure=(
            "A 9 to 12 residue loop at the second position and no Loop 1. A "
            "Parabacteroides merdae enzyme is the type example."
        ),
        substrates=(
            "A minority class in the published surveys, at roughly 4 to 7 per cent of gut "
            "GUS sequences. Substrate behaviour is less characterised than for Loop 1."
        ),
        does_not_mean=(
            "Being less characterised is a gap in the literature, not a finding about "
            "this sample."
        ),
    ),
    GusClass(
        code="mL1,2", label="mini-Loop 1,2",
        structure="Short loops at both positions, and neither at full length.",
        substrates=(
            "Described in the structural survey as a distinct category. Individual "
            "kinetics have not been reported for most members."
        ),
        does_not_mean="Carrying both mini-loops is not the same as carrying both full loops.",
    ),
    GusClass(
        code="NL", label="No Loop",
        structure=(
            "No loop at either active-site position. This is the largest class in the "
            "human gut surveys, and it is the class that contains every FMN-binding GUS "
            "identified so far, including the enzyme from Faecalibacterium prausnitzii "
            "L2-6 whose structure was solved with FMN bound 30 angstroms from the active "
            "site, and the confirmed binders from Ruminococcus gnavus and Roseburia "
            "hominis."
        ),
        substrates=(
            "Without the active-site loop these enzymes are generally slower on the small "
            "glucuronides that Loop 1 enzymes handle. The FMN-binding subset is defined by "
            "a surface site that affects the Michaelis constant but is not required for "
            "activity."
        ),
        does_not_mean=(
            "This classification does not resolve which No Loop enzymes are FMN binders. "
            "That was determined by structure and binding measurement on purified protein, "
            "not by loop length, so no FMN-binding count is reported here."
        ),
    ),
    GusClass(
        code="NC", label="No coverage",
        structure=(
            "The reference aligned to the E. coli enzyme, but the alignment did not span "
            "both loop positions, so neither loop length could be measured."
        ),
        substrates="Not classifiable, and therefore not interpreted.",
        does_not_mean=(
            "This is a limit of the alignment, not a property of the enzyme. These "
            "fragments are counted in the total and shown here so that the classified "
            "percentages have an honest denominator."
        ),
    ),
)

BY_CODE: Final[Mapping[str, GusClass]] = {c.code: c for c in CLASSES}


# --------------------------------------------------------------------------- #
# classification
# --------------------------------------------------------------------------- #


def insertion_length(
    query_gapped: str, subject_gapped: str, subject_start: int, window: tuple[int, int],
    *, slack: int = 3,
) -> int | None:
    """Residues the query *inserts* relative to the anchor at a window.

    The right measurement where the anchor has no loop of its own, which is
    the case at the second position: the published definition speaks of "a
    >=12-residue insertion region 416-419", and E. coli GUS has four ordinary
    residues there. Counting occupancy instead would score every reference
    that merely aligns normally, and the count would then rise with the width
    of the window rather than with anything about the enzyme.

    ``slack`` widens the span a little, because an alignment is free to place
    an insertion at either end of a run of equivalent positions and the
    choice is arbitrary. Insertions are extra residues by construction, so
    unlike occupancy this does not inflate with the width.
    """
    if len(query_gapped) != len(subject_gapped):
        return None
    lo, hi = window
    lo, hi = lo - slack, hi + slack
    s_pos = subject_start
    covered_lo = covered_hi = False
    count = 0
    inside = False
    for q_char, s_char in zip(query_gapped, subject_gapped, strict=True):
        if s_char != GAP:
            if s_pos <= lo:
                covered_lo = True
            if s_pos >= hi:
                covered_hi = True
            inside = lo <= s_pos <= hi
            s_pos += 1
        elif inside and q_char != GAP:
            count += 1
    if not (covered_lo and covered_hi):
        return None
    return count


def occupancy_length(
    query_gapped: str, subject_gapped: str, subject_start: int, window: tuple[int, int]
) -> int | None:
    """Residues of the query occupying an anchor window.

    The right measurement where the anchor carries the loop itself, which is
    the case at the first position: E. coli GUS has a 25-residue loop at
    356-380, so an enzyme without one aligns gaps across the window and
    scores near zero. Counts query residues aligned within the window and any
    the query inserts into it. ``None`` when the alignment does not cover the
    window, which is reported rather than treated as a loop of length zero —
    an enzyme whose alignment stops short is not an enzyme without a loop.
    """
    if len(query_gapped) != len(subject_gapped):
        return None
    lo, hi = window
    s_pos = subject_start
    covered_lo = covered_hi = False
    count = 0
    inside = False
    for q_char, s_char in zip(query_gapped, subject_gapped, strict=True):
        if s_char != GAP:
            if s_pos == lo:
                covered_lo = True
            if s_pos == hi:
                covered_hi = True
            inside = lo <= s_pos <= hi
            if inside and q_char != GAP:
                count += 1
            s_pos += 1
        elif inside and q_char != GAP:
            # An insertion relative to the anchor, sitting inside the window.
            count += 1
    if not (covered_lo and covered_hi):
        return None
    return count


def classify(loop1: int | None, loop2: int | None) -> str:
    """Assign a structural class from the two measured loop lengths."""
    if loop1 is None or loop2 is None:
        return "NC"
    if loop1 >= L1_FULL_MIN:
        return "L1"
    mini1 = L1_MINI_MIN <= loop1 < L1_FULL_MIN
    # The published ranges touch at 12 residues, which is both the floor of
    # Loop 2 and the ceiling of mini-Loop 2. Resolved in favour of the full
    # loop, so the class described as carrying extra binding domains is not
    # quietly shrunk by a boundary convention.
    full2 = loop2 >= L2_FULL_MIN
    mini2 = L2_MINI_MIN <= loop2 < L2_FULL_MIN
    if mini1 and (mini2 or full2):
        return "mL1,2"
    if mini1:
        return "mL1"
    if full2:
        return "L2"
    if mini2:
        return "mL2"
    return "NL"


def read_fasta(path: Path) -> Iterable[tuple[str, str]]:
    name: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].strip(), []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def gus_references(combined_faa: Path) -> dict[str, str]:
    """Accession to sequence, for the GUS entry only."""
    marker = f"~{GUS_ENTRY_ID}~"
    out: dict[str, str] = {}
    for header, sequence in read_fasta(combined_faa):
        if marker in header:
            out[header.rsplit("~", 1)[-1]] = sequence
    return out


def classify_references(
    combined_faa: Path, *, cache_dir: Path, diamond: str = "diamond", threads: int = 4,
) -> dict[str, Any]:
    """Classify every GUS reference, caching on the reference set itself."""
    refs = gus_references(combined_faa)
    if not refs:
        raise GusClassError(
            f"no references for {GUS_ENTRY_KEY} in {combined_faa}; the beta-glucuronidase "
            "panel is not in this database"
        )
    if ANCHOR_ACCESSION not in refs:
        raise GusClassError(
            f"the anchor {ANCHOR_ACCESSION} (E. coli GUS) is not among the "
            f"{len(refs)} references. Every loop position in the literature is numbered "
            "against it, so without it nothing here can be classified."
        )
    digest = hashlib.sha256(
        "|".join(f"{a}:{len(s)}" for a, s in sorted(refs.items())).encode()
    ).hexdigest()[:16]
    cache = cache_dir / f"gus_classes_{digest}.json"
    if cache.is_file():
        return json.loads(cache.read_text(encoding="utf-8"))

    exe = shutil.which(diamond)
    if exe is None:
        raise GusClassError(f"{diamond!r} not found on PATH")
    cache_dir.mkdir(parents=True, exist_ok=True)
    work = cache_dir / f"gus_work_{digest}"
    work.mkdir(parents=True, exist_ok=True)
    anchor_faa, refs_faa = work / "anchor.faa", work / "refs.faa"
    anchor_faa.write_text(
        f">{ANCHOR_ACCESSION}\n{refs[ANCHOR_ACCESSION]}\n", encoding="utf-8")
    refs_faa.write_text(
        "".join(f">{a}\n{s}\n" for a, s in sorted(refs.items())), encoding="utf-8")
    db, out_tsv = work / "anchor", work / "hits.tsv"
    _run([exe, "makedb", "--in", str(anchor_faa), "--db", str(db),
          "--threads", str(threads), "--quiet"])
    _run([exe, "blastp", "--db", str(db), "--query", str(refs_faa), "--out", str(out_tsv),
          "--outfmt", "6", "qseqid", "sstart", "send", "qseq_gapped", "sseq_gapped",
          "--max-target-seqs", "1", "--index-chunks", "1", "--very-sensitive",
          "--evalue", "1e-3", "--threads", str(threads), "--quiet"])

    assignment: dict[str, str] = {}
    lengths: dict[str, list[int | None]] = {}
    with out_tsv.open("r", encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            accession = parts[0]
            if accession in assignment:
                continue  # first (best) hit only
            try:
                s_start, s_end = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            lo = min(s_start, s_end)
            loop1 = occupancy_length(parts[3], parts[4], lo, L1_WINDOW)
            loop2 = insertion_length(parts[3], parts[4], lo, L2_WINDOW)
            assignment[accession] = classify(loop1, loop2)
            lengths[accession] = [loop1, loop2]

    for accession in refs:
        assignment.setdefault(accession, "NC")

    # ---- self-test ------------------------------------------------------- #
    # E. coli GUS is the enzyme the published definition is written against
    # and is a Loop 1 enzyme with a 25-residue loop. If the alignment walk is
    # off, this is where it shows, and every other assignment would be wrong
    # by the same amount.
    self_call = assignment.get(ANCHOR_ACCESSION)
    self_len = (lengths.get(ANCHOR_ACCESSION) or [None, None])[0]
    if self_call != "L1":
        raise GusClassError(
            f"the anchor {ANCHOR_ACCESSION} classified itself as {self_call!r} with a "
            f"Loop 1 length of {self_len}. E. coli GUS is the type example of Loop 1 "
            "with a 25-residue loop, so the alignment walk is wrong and no assignment "
            "here can be trusted."
        )

    counts: dict[str, int] = {c.code: 0 for c in CLASSES}
    for code in assignment.values():
        counts[code] = counts.get(code, 0) + 1
    payload = {
        "method": METHOD_GUS_CLASSES,
        "anchor": ANCHOR_ACCESSION,
        "anchor_loop1_residues": self_len,
        "reference_digest": digest,
        "n_references": len(refs),
        "reference_counts": counts,
        "assignment": assignment,
        "loop_lengths": dict(lengths),
    }
    cache.write_text(json.dumps(payload) + "\n", encoding="utf-8")
    for leftover in work.glob("*"):
        leftover.unlink(missing_ok=True)
    work.rmdir()
    return payload


def _run(cmd: Sequence[str]) -> None:
    try:
        subprocess.run(list(cmd), check=True, capture_output=True, text=True, timeout=1800)  # noqa: S603
    except subprocess.CalledProcessError as exc:
        raise GusClassError(
            f"command failed: {' '.join(cmd)}\n  {(exc.stderr or '').strip()[:400]}"
        ) from exc
    except subprocess.TimeoutExpired as exc:
        raise GusClassError(f"command timed out: {' '.join(cmd)} ({exc})") from exc


# --------------------------------------------------------------------------- #
# the reading
# --------------------------------------------------------------------------- #


#: The enzymes whose structural class is published, with the organism string
#: to find them by. A classifier that disagrees with the papers' own type
#: examples is not measuring what the papers measured, whatever else it is
#: doing, so this is checked before any breakdown is shown.
TYPE_EXAMPLES: Final[tuple[tuple[str, str], ...]] = (
    ("Escherichia coli", "L1"),
    ("Clostridium perfringens", "L1"),
    ("Streptococcus agalactiae", "L1"),
    ("Bacteroides fragilis", "mL1"),
    ("Bacteroides uniformis", "L2"),
    ("Parabacteroides merdae", "mL2"),
)

#: How many type examples must agree before a breakdown is reported. All of
#: them: these are six enzymes with solved structures, named in the source as
#: the definition of each class.
REQUIRED_AGREEMENT: Final = len(TYPE_EXAMPLES)


def validate(
    classification: Mapping[str, Any], organisms: Mapping[str, str]
) -> dict[str, Any]:
    """Check the assignment against the published type examples.

    Returns the evidence either way. A reference set contains many proteins
    per species, so a type example counts as agreeing when *any* of that
    species' references carries the published class — a deliberately generous
    test, because a species may genuinely encode several glucuronidases in
    different classes.
    """
    assignment = classification.get("assignment") or {}
    checks: list[dict[str, Any]] = []
    for needle, expected in TYPE_EXAMPLES:
        found = sorted({
            assignment[accession]
            for accession, organism in organisms.items()
            if needle in organism and accession in assignment
        })
        checks.append({
            "organism": needle, "expected": expected, "observed": found,
            "present": bool(found), "agrees": expected in found,
        })
    n_agree = sum(1 for c in checks if c["agrees"])
    return {
        "checks": checks,
        "n_agree": n_agree,
        "n_examples": len(checks),
        "passed": n_agree >= REQUIRED_AGREEMENT,
    }


#: Why the breakdown is withheld when validation fails. Written once, here,
#: so the report and the coverage ledger say the same thing.
WITHHELD_REASON: Final = (
    "The structural class of a beta-glucuronidase is defined by the length of a loop "
    "at one of two active-site positions, measured in a curated alignment guided by "
    "solved structures. Reproducing that from a pairwise alignment to the E. coli "
    "enzyme works for enzymes close to it and fails for the divergent gut ones: in "
    "testing, the Bacteroides fragilis enzyme came out as Loop 1 where the published "
    "structure is mini-Loop 1, and the Bacteroides uniformis enzyme came out as No "
    "Loop where it is the published example of Loop 2. Those are the gut organisms "
    "the breakdown would exist to describe, so no breakdown is shown. Your total "
    "beta-glucuronidase reading is unaffected and is reported as usual."
)

WHAT_WOULD_UNLOCK: Final = (
    "The curated loop alignment published as Data S3 of the source, or a profile "
    "model built from the solved structures, in place of pairwise alignment to a "
    "single anchor. The per-reference counts this needs are already recorded, so the "
    "breakdown appears as soon as the assignment agrees with the published examples."
)


def substrate_table(
    results: Mapping[str, Any],
    classification: Mapping[str, Any],
    *,
    validation: Mapping[str, Any] | None = None,
) -> dict[str, Any] | None:
    """Fragments by structural class, for this sample."""
    entry = _gus_entry(results)
    if entry is None:
        return None
    if validation is not None and not validation.get("passed"):
        return {
            "method": METHOD_GUS_CLASSES,
            "available": False,
            "reason": WITHHELD_REASON,
            "what_would_unlock_it": WHAT_WOULD_UNLOCK,
            "validation": dict(validation),
            "classes": [],
        }
    pairs = entry.get("accession_fragments")
    if not pairs:
        return {
            "method": METHOD_GUS_CLASSES,
            "available": False,
            "reason": (
                "this sample was screened before per-reference counts were retained, so "
                "its beta-glucuronidase fragments cannot be attributed to a structural "
                "class. Re-running the sample produces the breakdown."
            ),
            "classes": [],
        }
    assignment = classification.get("assignment") or {}
    ref_counts = classification.get("reference_counts") or {}
    by_class: dict[str, int] = {c.code: 0 for c in CLASSES}
    unknown = 0
    for accession, count in pairs:
        code = assignment.get(str(accession))
        if code is None:
            unknown += int(count)
            continue
        by_class[code] = by_class.get(code, 0) + int(count)
    total = sum(by_class.values()) + unknown
    rows = []
    for spec in CLASSES:
        n = by_class.get(spec.code, 0)
        rows.append({
            **spec.to_json(),
            "fragments": n,
            "percent_of_gus": round(100.0 * n / total, 1) if total else None,
            "reference_count": int(ref_counts.get(spec.code, 0)),
        })
    return {
        "method": METHOD_GUS_CLASSES,
        "available": True,
        "anchor": classification.get("anchor"),
        "total_fragments": total,
        "unassigned_fragments": unknown,
        "n_references": classification.get("n_references"),
        "classes": rows,
        "dominant": max(rows, key=lambda r: r["fragments"])["code"] if total else None,
        "caveats": [
            "The structural class is a property of each reference enzyme, established by "
            "solved structures and by loop length measured against the E. coli enzyme. "
            "It is not read off the sequencing read, which is far too short to see a "
            "whole protein.",
            "One species commonly encodes several glucuronidases in different classes, "
            "which is why this is attributed per reference protein rather than per "
            "organism.",
            "A class tells you which substrates an enzyme is plausibly fast on. It is not "
            "an activity measurement and cannot estimate exposure to any drug or hormone. "
            "Measured faecal beta-glucuronidase activity is a different assay, and can be "
            "entered in the lab results section if you have it.",
        ],
        "source": (
            "Pollet R.M. et al., \"An Atlas of beta-Glucuronidases in the Human "
            "Intestinal Microbiome\", Structure 25:967-977 (2017), "
            "doi:10.1016/j.str.2017.05.003; Little M.S. et al., \"Discovery and "
            "characterization of FMN-binding beta-glucuronidases in the human gut "
            "microbiome\", J Mol Biol 430:1-12 (2018), doi:10.1016/j.jmb.2018.12.013."
        ),
    }


def reference_organisms(index_json: Path) -> dict[str, str]:
    """Accession to organism, for the GUS entry, from the built database."""
    payload = json.loads(index_json.read_text(encoding="utf-8"))
    columns = payload.get("sequence_columns") or []
    try:
        i_key = columns.index("entry_key")
        i_acc = columns.index("accession")
        i_org = columns.index("organism")
    except ValueError:
        return {}
    return {
        row[i_acc]: row[i_org]
        for row in payload.get("sequences") or []
        if row[i_key] == GUS_ENTRY_KEY
    }


def analyse(
    results: Mapping[str, Any], *, refs_root: Path, diamond: str = "diamond",
) -> dict[str, Any] | None:
    """The A08 beta-glucuronidase substrate view for one sample.

    Returns ``None`` when the sample carries no beta-glucuronidase panel at
    all, and otherwise a view that either holds the breakdown or states, in
    the same shape, why it is withheld.
    """
    if _gus_entry(results) is None:
        return None
    fingerprint = str((results.get("run") or {}).get("reference fingerprint") or "")
    db_dir = refs_root / "db" / fingerprint
    combined, index_json = db_dir / "combined.faa", db_dir / "index.json"
    if not fingerprint or not combined.is_file() or not index_json.is_file():
        return {
            "method": METHOD_GUS_CLASSES,
            "available": False,
            "reason": (
                "the reference database this sample was screened against is no longer "
                "on disk, so its beta-glucuronidase references cannot be classified."
            ),
            "what_would_unlock_it": "rebuilding the database and re-running the sample",
            "classes": [],
        }
    classification = classify_references(
        combined, cache_dir=db_dir, diamond=diamond)
    checked = validate(classification, reference_organisms(index_json))
    return substrate_table(results, classification, validation=checked)


def _gus_entry(results: Mapping[str, Any]) -> Mapping[str, Any] | None:
    for panel in results.get("panels") or []:
        if not isinstance(panel, dict) or panel.get("name") != "bglucuronidase":
            continue
        for gene in panel.get("genes") or []:
            if isinstance(gene, dict) and gene.get("entry_id") == GUS_ENTRY_ID:
                return gene
    return None
