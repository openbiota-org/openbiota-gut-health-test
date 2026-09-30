"""The experimental P9 protein — spec 0.8.3 §9.1, source C02.

P9 is a single *Akkermansia muciniphila* protein with experimental evidence
for an effect on GLP-1 secretion. §9.1 requires it to be detected from the
exact sequence that was tested, and is emphatic about why: a species
abundance is not a substitute for a protein, and a divergent homolog of a
protease family is not the protein either.

So the detector is the tested sequence and nothing else. Everything in the
same family competes against it rather than counting towards it, which is
what "competitive mapping against related S41 proteases" means: a read from
some other carboxyl-terminal protease lands on the family decoy, not on P9.

Two identifiers are easy to confuse and the specification names the trap
explicitly. `Amuc_1631` is P9, 748 residues. `Amuc_1831` is a 98-residue
four-helix bundle protein and is not an accepted alias for anything here.
The guard against that mistake is a hash: both the protein and the coding
sequence are pinned by SHA-256, verified against the live records, and a
mismatch is an error rather than a warning.

What the reading says, and only this: how much read evidence maps to the
experimental sequence, over how much of its length, and how that compares
with the family around it. Not protein secretion, which needs proteomics.
Not a GLP-1 concentration, which needs an assay on blood with its timing
recorded. Not an equivalence to GLP-1 medication, which the specification
also prohibits and which would be absurd from a stool gene count.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Final

METHOD_P9: Final = "ext083.host_signaling.p9_sequence_potential/1.0"
METRIC_ID: Final = "host_signaling.p9_sequence_potential"

#: The panel entry that carries the reading, and the family it competes with.
PANEL: Final = "p9"
TARGET_ENTRY: Final = "P9"
FAMILY_DECOY: Final = "S41A"
WRONG_ID_DECOY: Final = "WRONGID"


class P9Error(ValueError):
    """An anchor that did not match what the specification pins."""


@dataclass(frozen=True)
class Anchor:
    """One pinned sequence record.

    ``sha256`` is over the uppercase, unwrapped sequence with no trailing
    newline; the coding sequence is hashed in its coding orientation, which
    for this locus means the reverse complement of the genome slice.
    """

    label: str
    accession: str
    url: str
    sha256: str
    length: int
    note: str = ""

    def verify(self, sequence: str) -> None:
        got = hashlib.sha256(sequence.upper().encode()).hexdigest()
        if len(sequence) != self.length:
            raise P9Error(
                f"{self.label} ({self.accession}): {len(sequence)} residues, expected "
                f"{self.length}. The record has changed or the wrong one was fetched."
            )
        if got != self.sha256:
            raise P9Error(
                f"{self.label} ({self.accession}): sha256 {got}, expected {self.sha256}. "
                "The sequence is not the one this detector was built against, so the "
                "detector cannot be trusted."
            )


PROTEIN: Final = Anchor(
    label="P9 protein",
    accession="B2UM07",
    url="https://rest.uniprot.org/uniprotkb/B2UM07.fasta",
    sha256="8618a2f5066940d6ae5bd5f99859f478a6dbc3036ec3b77eb59cc3412160f6a2",
    length=748,
    note=(
        "Carboxyl-terminal protease, locus Amuc_1631, GenBank ACD05451.1. Belongs to "
        "the peptidase S41A family, which is why the family is a competitor here rather "
        "than corroboration."
    ),
)

#: The genome record and the coordinates the coding sequence is taken from.
GENOME_ACCESSION: Final = "CP001071.1"
GENOME_URL: Final = "https://www.ebi.ac.uk/ena/browser/api/fasta/CP001071.1"
CDS_SPAN: Final[tuple[int, int]] = (1_965_361, 1_967_607)
CDS_IS_COMPLEMENT: Final = True

CDS: Final = Anchor(
    label="P9 coding sequence",
    accession=f"{GENOME_ACCESSION}:complement({CDS_SPAN[0]}..{CDS_SPAN[1]})",
    url=GENOME_URL,
    sha256="fb313b8d177bed36da19ae85884357114f00e29cf47387018b3d39fc55f95bdc",
    length=2_247,
    note="2,247 nucleotides including the stop codon, in coding orientation.",
)

#: The identifier the specification names as a wrong-identifier control. It
#: is a different protein of a different length, and is not an alias.
WRONG_IDENTIFIER: Final[Mapping[str, Any]] = {
    "locus": "Amuc_1831",
    "genbank": "ACD05649.1",
    "uniprot": "B2UN38",
    "length": 98,
    "protein": "Four helix bundle protein",
    "why_it_is_here": (
        "Published discussion of P9 has cited this locus by mistake. It is 98 residues "
        "against P9's 748 and belongs to no protease family, so a detector that accepted "
        "it would be measuring something unrelated. The identity fixture requires it to "
        "fail."
    ),
}


def reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def extract_cds(genome: str) -> str:
    """The coding sequence, from the genome, in coding orientation."""
    start, end = CDS_SPAN
    if len(genome) < end:
        raise P9Error(
            f"{GENOME_ACCESSION}: {len(genome):,} nucleotides, too short to hold "
            f"{start:,}..{end:,}. The wrong record was fetched."
        )
    piece = genome[start - 1:end]
    return reverse_complement(piece) if CDS_IS_COMPLEMENT else piece.upper()


def verify_anchors(*, cache_dir: Path, fetch: Any = None) -> dict[str, Any]:
    """Check both anchors against the pinned hashes, and record the outcome.

    Cached on success so a report does not depend on two networks being up,
    and so the recorded provenance is the thing that was actually verified
    rather than a claim about it.
    """
    cache = cache_dir / "p9_anchors.json"
    if cache.is_file():
        return json.loads(cache.read_text(encoding="utf-8"))
    if fetch is None:
        import urllib.request

        def fetch(url: str) -> str:  # noqa: ANN001 - local default
            request = urllib.request.Request(url, headers={"User-Agent": "openbiota"})
            with urllib.request.urlopen(request, timeout=300) as handle:  # noqa: S310
                return handle.read().decode()

    protein = "".join(fetch(PROTEIN.url).splitlines()[1:]).upper()
    PROTEIN.verify(protein)

    genome = "".join(
        line.strip() for line in fetch(GENOME_URL).splitlines()
        if not line.startswith(">")
    ).upper()
    cds = extract_cds(genome)
    CDS.verify(cds)

    payload = {
        "method": METHOD_P9,
        "protein": {
            "accession": PROTEIN.accession, "length": len(protein),
            "sha256": PROTEIN.sha256, "verified": True, "note": PROTEIN.note,
        },
        "cds": {
            "accession": CDS.accession, "length": len(cds),
            "sha256": CDS.sha256, "verified": True,
            "orientation": "coding (reverse complement of the genome slice)",
            "note": CDS.note,
        },
        "genome": {"accession": GENOME_ACCESSION, "length": len(genome)},
        "wrong_identifier": dict(WRONG_IDENTIFIER),
    }
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache.write_text(json.dumps(payload, indent=1) + "\n", encoding="utf-8")
    return payload


# --------------------------------------------------------------------------- #
# the reading
# --------------------------------------------------------------------------- #


def reading(
    results: Mapping[str, Any], *, anchors: Mapping[str, Any] | None = None
) -> dict[str, Any]:
    """`host_signaling.p9_sequence_potential` for one sample."""
    panel = next(
        (p for p in results.get("panels") or []
         if isinstance(p, dict) and p.get("name") == PANEL),
        None,
    )
    if panel is None:
        return {
            "method": METHOD_P9,
            "metric_id": METRIC_ID,
            "state": "not_assayed",
            "reason": (
                "the P9 sequence panel is not in the database this sample was screened "
                "against, so the experimental protein was not searched for"
            ),
        }

    target = next(
        (g for g in panel.get("genes") or [] if g.get("entry_id") == TARGET_ENTRY), None)
    decoys = {d.get("entry_id"): d for d in panel.get("decoys") or []}
    family = decoys.get(FAMILY_DECOY) or {}
    wrong = decoys.get(WRONG_ID_DECOY) or {}

    fragments = int((target or {}).get("fragments") or 0)
    family_fragments = int(family.get("fragments") or 0)
    identity = (target or {}).get("median_identity")
    state = "present" if fragments > 0 else "absent"

    total = fragments + family_fragments
    return {
        "method": METHOD_P9,
        "metric_id": METRIC_ID,
        "state": state,
        "fragments": fragments,
        "median_identity": identity,
        "copies_per_100_genomes": (target or {}).get("copies_per_100_genomes"),
        "reference_count": (target or {}).get("reference_count"),
        # §9.1: competitive mapping against the related S41 proteases, so a
        # read from the family cannot be read as the experimental protein.
        "family_competition": {
            "entry": FAMILY_DECOY,
            "fragments": family_fragments,
            "share_of_family_and_target": (
                round(fragments / total, 4) if total else None
            ),
            "meaning": (
                "Fragments on the wider carboxyl-terminal protease family, which this "
                "protein belongs to. They are shown beside the target and never added "
                "to it: a homolog of a protease family is not the experimental protein."
            ),
        },
        "wrong_identifier_control": {
            "entry": WRONG_ID_DECOY,
            "fragments": int(wrong.get("fragments") or 0),
            "meaning": WRONG_IDENTIFIER["why_it_is_here"],
        },
        "anchors": dict(anchors or {}),
        "establishes": (
            "read evidence mapping to the exact protein sequence that was tested "
            "experimentally, at the panel's identity floor"
        ),
        "does_not_establish": (
            "that the protein is being made, that it is being secreted, or anything "
            "about your GLP-1. Secretion needs proteomics and a hormone level needs a "
            "blood assay with its timing recorded. This is also not comparable to a "
            "GLP-1 medication, and nothing here should be read as though it were."
        ),
        "limitations": (
            "One reference sequence means low sensitivity by design: a diverged P9 in "
            "another strain will not reach the identity floor and is reported as "
            "unresolved rather than transferred onto this annotation.",
            "The evidence for an effect is preclinical, from mice and from cell work.",
        ),
        "source_ids": ("C02", "F43"),
    }
