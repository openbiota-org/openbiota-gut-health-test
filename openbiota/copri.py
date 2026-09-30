"""The *Segatella copri* complex, resolved to clade level.

*Prevotella copri* is the organism most often named in the rheumatoid-arthritis
microbiome literature, and for years it was also the least interpretable thing
a species-level report could print. Tett et al. showed in 2019 that "P. copri"
is not one species but four deeply divergent clades; Manghi et al. extended
that to 13 in 2023, and the genus was renamed *Segatella* along the way. A
single "P. copri: 2.1%" line silently averaged organisms 13–21% divergent from
one another — further apart than many separate genera.

Two things follow, and this module exists to make both of them concrete.

**The clades are measurable, and we already measure them.** Each clade is a
distinct SGB in the MetaPhlAn 4 catalogue the extended lane already runs. So a
report does not have to say "strain resolution unavailable": it can name which
clades are present, at what abundance, and — when none is — say so definitely
across every clade assessed. A definite negative is worth more than an
ambiguity.

**No clade has an established disease direction.** Manghi et al. tested all 13
against 1,635 cases and 1,854 controls over 22 studies and 11 conditions,
including a rheumatoid-arthritis cohort, and found none of the complex's
species associated with any condition. So this module reports composition and
nothing else. It scores no disease resemblance, because there is no published
clade-level direction to score against, and inventing one would be worse than
saying nothing.

Sources:

* Tett A, et al. The Prevotella copri complex comprises four distinct clades
  underrepresented in Westernized populations. Cell Host Microbe
  2019;26(5):666-679. doi:10.1016/j.chom.2019.08.018
* Manghi P, et al. Extension of the Segatella copri complex to 13 species with
  distinct large extrachromosomal elements and associations with host
  conditions. Cell Host Microbe 2023;31(11):1856-1872.
  doi:10.1016/j.chom.2023.09.013
* Hitch TCA, et al. A taxonomic note on the genus Prevotella. Syst Appl
  Microbiol 2022;45:126354. doi:10.1016/j.syapm.2022.126354
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Final


@dataclass(frozen=True, slots=True)
class Clade:
    """One species-level clade of the complex."""

    clade: str
    species: str
    #: SGB identifiers seen for this clade. Matched alongside the species name
    #: because MetaPhlAn renames species between database releases while SGB
    #: numbers stay put.
    sgbs: tuple[str, ...]
    host: str
    note: str = ""


#: Clade letters from Tett 2019 (A–D) and Manghi 2023 (E–M), with the species
#: names Manghi proposed and the SGBs they carry in the CHOCOPhlAn SGB
#: catalogue. Clades J–M were recovered only from non-human primates and are
#: listed so that "assessed" means the whole complex, not just the human part.
CLADES: Final[tuple[Clade, ...]] = (
    Clade("A", "Segatella_copri", ("SGB1626",), "human",
          "The clade that carries the original P. copri name and most "
          "Westernized-population genomes."),
    Clade("B", "Segatella_brunsvicensis", ("SGB1613",), "human"),
    Clade("C", "Segatella_sinensis", ("SGB1644", "SGB1638"), "human",
          "SGB1638 is spelt Sgatella_sinensis in some catalogue releases."),
    Clade("D", "Segatella_brasiliensis", ("SGB1635", "SGB1636"), "human"),
    Clade("E", "Segatella_hominis", ("SGB1615", "SGB1617"), "human",
          "Formerly Prevotella hominis."),
    Clade("F", "Segatella_sanihominis", ("SGB1614",), "human"),
    Clade("G", "Segatella_sinica", ("SGB1653",), "human"),
    Clade("H", "Candidatus_Segatella_caccae", ("SGB1612",), "human"),
    Clade("I", "Candidatus_Segatella_intestinihominis", ("SGB1623",), "human"),
    Clade("J", "Candidatus_Segatella_violae", ("SGB20122",), "non_human_primate"),
    Clade("K", "Candidatus_Segatella_albertsiae", ("SGB20121",), "non_human_primate"),
    Clade("L", "Candidatus_Segatella_mututuai", (), "non_human_primate"),
    Clade("M", "Candidatus_Segatella_papionis", (), "non_human_primate"),
)

#: What the literature does and does not support at clade level. Printed with
#: the composition so a reader never has to guess whether a present clade is
#: the "bad" one — the answer is that no such clade has been established.
DISEASE_DIRECTION: Final = (
    "No clade of this complex has an established association with any disease, "
    "including rheumatoid arthritis. Manghi et al. tested all 13 clades against "
    "1,635 disease cases and 1,854 controls across 22 studies and 11 conditions, "
    "a set that included a rheumatoid arthritis cohort, and found none of the "
    "complex's species associated with any health condition. Earlier "
    "rheumatoid-arthritis reports named the species as a whole, not a clade, and "
    "those reports have replicated poorly: in an anti-CCP-positive at-risk cohort "
    "three Prevotellaceae strains were enriched and five depleted in the people "
    "who progressed. So this section reports what is present and makes no disease "
    "claim from it."
)

SOURCES: Final = (
    "Tett 2019 Cell Host Microbe doi:10.1016/j.chom.2019.08.018",
    "Manghi 2023 Cell Host Microbe doi:10.1016/j.chom.2023.09.013",
    "Hitch 2022 Syst Appl Microbiol doi:10.1016/j.syapm.2022.126354",
)


def _norm(text: str) -> str:
    return str(text or "").strip().lower().replace(" ", "_")


def resolve(sgbs: Iterable[Mapping[str, Any]] | None) -> dict[str, Any]:
    """Clade composition of the complex from the extended (SGB) lane.

    `sgbs` is `extended_catalogue["sgbs"]`: one row per detected SGB with
    `sgb`, `species` and `percent`. When the lane did not run, the result says
    so rather than reporting an absence nobody measured.
    """
    if sgbs is None:
        return {
            "status": "not_assessed",
            "reason": (
                "the extended SGB catalogue lane did not run, and the species-level lane "
                "cannot separate the clades of this complex"
            ),
            "clades_assessed": 0,
            "clades_present": [],
            "total_percent": None,
            "disease_direction": DISEASE_DIRECTION,
            "sources": list(SOURCES),
        }

    rows = list(sgbs)
    by_sgb: dict[str, Mapping[str, Any]] = {}
    by_species: dict[str, Mapping[str, Any]] = {}
    for row in rows:
        if sgb := str(row.get("sgb") or "").strip():
            by_sgb[sgb] = row
        if species := _norm(row.get("species")):
            by_species[species] = row

    present: list[dict[str, Any]] = []
    total = 0.0
    for clade in CLADES:
        hit = next((by_sgb[s] for s in clade.sgbs if s in by_sgb), None)
        if hit is None:
            hit = by_species.get(_norm(clade.species))
        if hit is None:
            continue
        percent = float(hit.get("percent") or 0.0)
        total += percent
        present.append(
            {
                "clade": clade.clade,
                "species": clade.species.replace("_", " "),
                "sgb": str(hit.get("sgb") or ""),
                "percent": round(percent, 4),
                "host_range": clade.host,
                "note": clade.note,
            }
        )

    human = [c for c in CLADES if c.host == "human"]
    return {
        "status": "resolved" if present else "absent",
        "clades_assessed": len(CLADES),
        "human_clades_assessed": len(human),
        "clades_present": present,
        "total_percent": round(total, 4),
        "plain": _plain(present, total, len(human), len(CLADES)),
        "disease_direction": DISEASE_DIRECTION,
        "what_this_replaces": (
            "A single species-level \u201cP. copri\u201d number, which averaged organisms "
            "13\u201321% divergent from one another and could not say which was present."
        ),
        "sources": list(SOURCES),
    }


def _plain(
    present: Sequence[Mapping[str, Any]], total: float, n_human: int, n_all: int
) -> str:
    if not present:
        return (
            f"No member of the Segatella (Prevotella) copri complex was detected. All "
            f"{n_all} clades were assessed, including the {n_human} found in humans, so this "
            "is a measured absence across the whole complex rather than an unresolved "
            "species-level reading."
        )
    named = ", ".join(
        f"clade {c['clade']} ({c['species']}, {c['percent']:.3f}%)" for c in present
    )
    return (
        f"The complex is present at {total:.3f}% of the community, resolved to "
        f"{len(present)} of {n_all} clades: {named}."
    )


def self_test() -> int:
    """Check the resolver against hand-built rows. Returns a failure count."""
    failures = 0

    absent = resolve([{"sgb": "SGB4285", "species": "Ruminococcus_bromii", "percent": 9.1}])
    if absent["status"] != "absent":
        failures += 1
    if absent["clades_assessed"] != len(CLADES):
        failures += 1
    if "measured absence" not in absent["plain"]:
        failures += 1

    # Matched on SGB id even when the species label has been renamed.
    renamed = resolve([{"sgb": "SGB1626", "species": "Prevotella_copri", "percent": 2.5}])
    if renamed["status"] != "resolved" or renamed["total_percent"] != 2.5:
        failures += 1
    if renamed["clades_present"][0]["clade"] != "A":
        failures += 1

    # And on species name when the SGB id is missing.
    by_name = resolve([{"sgb": "", "species": "Segatella hominis", "percent": 0.4}])
    if by_name["status"] != "resolved" or by_name["clades_present"][0]["clade"] != "E":
        failures += 1

    if resolve(None)["status"] != "not_assessed":
        failures += 1

    return failures


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(self_test())
