"""Reconcile the pathogen screen's bacterial calls with the organism inventory.

The screen aligns reads to its target genomes and decoys. A fragment it
counts as "unique" is unique among *those* references, not among the
organisms in the sample: a relative the bundle does not carry (*Dorea
hominis*, *Mediterraneibacter lactaris*, the Blautias) has nowhere else to
put its conserved reads, so they land on the nearest target and read as a
signal. One sample showed 19,156 fragments on *M. gnavus* - "0.34% of
analysed DNA" - while nine independent methods, mapping the same reads
competitively against every relative they had found, saw no *M. gnavus* at
all. The coverage gave it away: 8% of the genome, piled into 106 conserved
regions at 98% identity, where a genuine population at that depth covers
four-fifths of its genome at 100%.

The organism inventory is the report's one account of which bacteria are
present. This pass makes the screen agree with it:

* a bacterial target the inventory **found** keeps its call and quotes the
  community share beside the screen's own figure;
* a target the inventory did **not** find, whose coverage has the shape of
  shared sequence (piled into conserved regions - low evenness - or below
  the identity a genuine population maps at) while relatives *were* found,
  is reported as **shared sequence from those relatives**, named, and no
  longer counts as a finding;
* a target the inventory did not find but whose coverage looks genuine is
  kept, and says that the screen's own alignment is its only evidence.

Nothing is deleted: every record keeps its fragments, regions, breadth and
identity, and gains the reconciliation fields below.
"""

from __future__ import annotations

import math
import re
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Final

from openbiota.expansion.confirm import MIN_EVENNESS, _evenness
from openbiota.inventory import Inventory, Organism, canonical, display_of

#: A screen count more than this many times the community share is padded
#: with relatives' shared sequence; the card says so and quotes the share.
SHARE_RATIO_NOTE: Final = 5.0
#: Fragment length assumed when the run's read length is unknown.
DEFAULT_FRAGMENT_BP: Final = 300
#: Sequence statuses the screen counts as a named finding.
_SUPPORTED: Final = frozenset({"supported_sequence", "marker_signal"})
_ATTENTION_CLASSES: Final = frozenset({
    "established_enteric", "rare_enteric", "toxin_or_pathotype_dependent", "extraintestinal_watch",
})


def _norm(name: str) -> str:
    return re.sub(r"[^a-z0-9]", "", name.lower())


def species_names(catalog_dir: Path | None = None) -> dict[str, str]:
    """target_id -> the species the target resolves to (the catalogue's
    ``resolution_name``: *Bacteroides fragilis* for the ETBF pathotype target,
    *Bacillus cereus* for the sensu-lato group)."""
    try:
        from openbiota.pathogens.catalog import load_catalog

        cat = load_catalog(catalog_dir or Path("pathogens/catalog"))
    except Exception:  # noqa: BLE001 - the catalogue is optional here; display names still work
        return {}
    return {tid: str(seed.resolution_name) for tid, seed in cat.by_id.items() if getattr(seed, "resolution_name", None)}


def _names_of(display_name: str) -> list[str]:
    """Every spelling the inventory might hold a target under."""
    base = re.sub(r"\s*\(.*?\)\s*", " ", display_name).strip()
    base = re.sub(r"^(enterotoxigenic|toxigenic|enteropathogenic|enteroaggregative|enteroinvasive|shiga toxin-producing)\s+",
                  "", base, flags=re.I)
    base = re.sub(r"\s+(sensu lato|sensu stricto|group|complex)$", "", base, flags=re.I).strip()
    under = base.replace(" ", "_")
    out = [under, canonical(under), base]
    try:
        from openbiota.expansion import names as _names

        bridged = _names.ncbi_species_to_r232().get(_names.mpa_style(canonical(under))) or ""
        if bridged:
            out.extend([bridged, bridged.replace(" ", "_")])
    except Exception:  # noqa: BLE001 - a lookup table; its absence must not cost the pass
        pass
    return list(dict.fromkeys(n for n in out if n))


def _lookup(inv: Inventory, display_name: str) -> Organism | None:
    for name in _names_of(display_name):
        o = inv.get(name)
        if o is not None:
            return o
    wanted = {_norm(n) for n in _names_of(display_name)}
    for o in inv.organisms:
        keys = {_norm(o.display), _norm(o.species), _norm(o.gtdb or ""), *(_norm(a) for a in o.aliases)}
        if keys & wanted:
            return o
    return None


def _found(o: Organism | None) -> bool:
    if o is None:
        return False
    return (o.in_primary and o.percent > 0) or (o.status == "supported" and (o.best_percent or 0.0) > 0)


def _rejected_match(rejected: list[Mapping[str, Any]], display_name: str) -> Mapping[str, Any] | None:
    wanted = {_norm(n) for n in _names_of(display_name)}
    for r in rejected:
        keys = {_norm(display_of(r)), _norm(str(r.get("species") or "")), _norm(str(r.get("gtdb") or "")),
                *(_norm(str(a)) for a in (r.get("aliases") or ()))}
        if keys & wanted:
            return r
    return None


def _genus_of(name: str) -> str:
    return name.replace("_", " ").split(" ")[0].split("_")[0]


def _relatives(inv: Inventory, display_name: str) -> list[Organism]:
    """Organisms found in the composition that share the target's genus - by
    GTDB placement where the bridge knows the name - or that an older
    catalogue listed under the target's name."""
    names = _names_of(display_name)
    genera = {_genus_of(n) for n in names}
    gtdb_genera = {n.split(" ")[0] for n in names if " " in n and n.split(" ")[0] != display_name.split(" ")[0]}
    wanted = {_norm(n) for n in names}
    out: list[Organism] = []
    for o in inv.organisms:
        if not (o.in_primary and o.percent > 0):
            continue
        listed = {_norm(x) for x in (o.formerly_listed_as or ())}
        genus_hit = ((o.gtdb_genus or (o.gtdb or "").split(" ")[0]) in gtdb_genera
                     or _genus_of(o.gtdb_genus or o.genus or o.display) in genera)
        if listed & wanted or genus_hit:
            out.append(o)
    # the population an older catalogue called by this very name first, then by share
    out.sort(key=lambda o: (not ({_norm(x) for x in (o.formerly_listed_as or ())} & wanted), -o.percent))
    return out


def evenness_of(record: Mapping[str, Any], fragment_bp: int = DEFAULT_FRAGMENT_BP) -> float | None:
    """Observed breadth over the breadth the same depth would give at random."""
    try:
        frags = int(record.get("unique_supporting_fragments") or 0)
        bases = float(record.get("informative_bases_covered") or 0.0)
        frac = float(record.get("informative_region_breadth_fraction") or record.get("reference_breadth_fraction") or 0.0)
    except (TypeError, ValueError):
        return None
    if frags <= 0 or bases <= 0 or frac <= 0:
        return None
    length = bases / frac
    depth = frags * fragment_bp / length
    breadth = float(record.get("reference_breadth_fraction") or frac)
    return _evenness(breadth, depth)


def _pct(v: float) -> str:
    return f"{v:.3f}%" if v < 1 else f"{v:.2f}%"


def reconcile(pathogens: dict[str, Any], inv: Inventory | None, *, read_length_bp: float | None = None,
              rejected: list[Mapping[str, Any]] | None = None,
              species_of: Mapping[str, str] | None = None) -> dict[str, Any]:
    """Mutate ``pathogens`` (the screen's JSON) so its bacterial calls agree
    with the inventory; return the reconciliation summary also written under
    ``pathogens["inventory_reconciliation"]``."""
    summary: dict[str, Any] = {
        "n_checked": 0, "n_found": 0, "n_shared_sequence": 0, "n_not_in_inventory": 0,
        "shared_sequence": [], "not_in_inventory": [],
        "rule": ("Bacterial calls with enough sequence to name the species are checked against the organism "
                 "inventory, which maps the same reads competitively against every relative it found. A species "
                 "the inventory did not find, whose coverage is piled into conserved regions (evenness below "
                 f"{MIN_EVENNESS}: a fraction of the genome at a depth that would have covered most of it), is "
                 "no longer counted as a finding. Absolute identity is not used: a genuine strain 1-2% divergent "
                 "from the reference maps at the same identity as a cross-mapped relative, so only the shape of "
                 "the coverage and the competitive evidence decide."),
    }
    if not isinstance(pathogens, dict) or inv is None:
        pathogens["inventory_reconciliation"] = {**summary, "status": "inventory unavailable"}
        return pathogens["inventory_reconciliation"]
    fragment_bp = int(round(2 * read_length_bp)) if read_length_bp else DEFAULT_FRAGMENT_BP
    rejected = list(rejected or [])
    species_of = dict(species_of) if species_of is not None else species_names()
    for rec in pathogens.get("results") or []:
        if not isinstance(rec, dict) or rec.get("group") != "bacteria":
            continue
        if rec.get("sequence_status") not in _SUPPORTED:
            continue
        if rec.get("species_resolution") not in (None, "resolved", "resolved_by_marker"):
            continue
        shown = str(rec.get("display_name") or rec.get("target_id") or "")
        name = species_of.get(str(rec.get("target_id") or ""), shown)
        summary["n_checked"] += 1
        rec["inventory_species"] = name
        o = _lookup(inv, name)
        identity = float(rec.get("median_alignment_identity") or 0.0)
        evenness = evenness_of(rec, fragment_bp)
        rec["inventory_evenness"] = None if evenness is None else round(evenness, 3)
        if _found(o):
            summary["n_found"] += 1
            rec["inventory_agreement"] = "found"
            rec["inventory_organism"] = o.display
            rec["inventory_share_percent"] = round(o.percent, 4) if o.in_primary else None
            rec["inventory_methods"] = len(o.methods or ())
            share = (f"In the community composition it is {_pct(o.percent)} (organism inventory"
                     f"{', seen by ' + str(len(o.methods)) + ' methods' if len(o.methods or ()) > 1 else ''})."
                     if o.in_primary and o.percent > 0 else
                     f"The organism inventory found it too ({len(o.methods or ())} methods).")
            fpm = rec.get("normalized_fragments_per_million")
            screen_pct = (float(fpm) / 10_000.0) if fpm is not None else None
            if (screen_pct is not None and o.in_primary and o.percent > 0
                    and screen_pct > SHARE_RATIO_NOTE * o.percent):
                share += (" The screen's own count is larger because it is not competitive: where relatives are "
                          "present their shared sequence lands on this target too, so the community share is the "
                          "amount.")
                rec.setdefault("reason_codes", [])
                if "count_includes_shared_sequence" not in rec["reason_codes"]:
                    rec["reason_codes"] = [*rec["reason_codes"], "count_includes_shared_sequence"]
            rec["plain_statement"] = (str(rec.get("plain_statement") or "").rstrip() + " " + share).strip()
            continue
        rej = _rejected_match(rejected, name)
        kin = _relatives(inv, name)
        # Only the coverage's own shape, or the competitive confirmation's
        # verdict, may withdraw a call. The inventory's silence alone must
        # not: its lanes have detection floors the screen does not, so a
        # genuine pathogen at 0.001% is invisible to them and must still be
        # reported.
        shared = rej is not None or (evenness is not None and evenness < MIN_EVENNESS)
        if shared:
            summary["n_shared_sequence"] += 1
            named: list[str] = []
            if rej is not None and rej.get("reads_belong_to"):
                named.append(str(rej["reads_belong_to"]))
            named.extend(f"{k.display} ({_pct(k.percent)})" for k in kin[:3] if k.display not in named)
            listed_as = [k for k in kin if {x.replace("_", " ").lower() for x in (k.formerly_listed_as or ())}
                         & {n.replace("_", " ").lower() for n in _names_of(name)}]
            why = []
            if evenness is not None:
                why.append(f"the fragments cover {float(rec.get('reference_breadth_fraction') or 0) * 100:.1f}% of the "
                           f"genome in {int(rec.get('informative_regions_supported') or 0):,} conserved regions "
                           f"(evenness {evenness:.2f}; a genuine population at this depth covers its genome evenly)")
            if identity:
                why.append(f"they map at {identity:.1%} identity")
            if rej is not None:
                why.append("the competitive whole-genome confirmation rejected the call for this sample")
            statement = (
                f"Sequence matched {shown}, but the nine-method organism inventory, which maps the same reads "
                f"competitively against every relative it found, did not find {name}. "
                + ("; ".join(why).capitalize() + ". " if why else "")
                + ("This is shared sequence from relatives that are present: " + ", ".join(named) + "."
                   if named else
                   "This is the shape of shared sequence from a relative rather than of the organism; no organism "
                   "in the inventory carries the name, so the relative is one the catalogues hold above species "
                   "level or not at all.")
                + (f" The population an older catalogue listed as {name} in this sample is "
                   f"{listed_as[0].display} ({_pct(listed_as[0].percent)})." if listed_as else "")
                + " It is not counted as a finding."
            )
            rec.update({
                "inventory_agreement": "shared_sequence_from_relatives",
                "inventory_relatives": named,
                "sequence_status": "ambiguous_signal",
                "display_status": "ambiguous_signal",
                "display_qualifier": "shared sequence from relatives found in this sample; species not found by the organism inventory",
                "species_resolution": "not_supported_by_inventory",
                "counts_as_pathogen": False,
                "report_tier": "shared_sequence",
                "plain_statement": statement,
                "reason_codes": [*dict.fromkeys([*(rec.get("reason_codes") or []), "shared_sequence_from_relatives"])],
            })
            summary["shared_sequence"].append({"target": shown, "species": name, "relatives": named,
                                               "identity": round(identity, 4), "evenness": rec["inventory_evenness"]})
            continue
        summary["n_not_in_inventory"] += 1
        rec["inventory_agreement"] = "not_in_inventory"
        rec["reason_codes"] = [*dict.fromkeys([*(rec.get("reason_codes") or []), "not_in_organism_inventory"])]
        rec["plain_statement"] = (str(rec.get("plain_statement") or "").rstrip()
                                  + " The nine-method organism inventory did not name this species - its reference "
                                    "catalogues have a detection floor this screen does not - and the coverage here "
                                    "has the shape of a genuine population, so the finding stands on this screen's "
                                    "own alignment.").strip()
        summary["not_in_inventory"].append(shown)
    _recount(pathogens)
    pathogens["inventory_reconciliation"] = {**summary, "status": "completed"}
    return pathogens["inventory_reconciliation"]


def _recount(pathogens: dict[str, Any]) -> None:
    """Recompute the screen's headline counts from the records as they now stand
    (mirrors `PathogenBranchResult`'s properties)."""
    recs = [r for r in (pathogens.get("results") or []) if isinstance(r, dict)]
    supported = [r for r in recs if r.get("sequence_status") in _SUPPORTED and r.get("display_status") != "technical_artifact"]
    counts = pathogens.get("counts") if isinstance(pathogens.get("counts"), dict) else {}
    counts.update({
        "pathogen_count": sum(1 for r in supported if r.get("interpretation_class") != "background_or_decoy"),
        "supported": len(supported),
        "attention": sum(1 for r in supported if r.get("interpretation_class") in _ATTENTION_CLASSES),
        "opportunists": sum(1 for r in supported if r.get("interpretation_class") == "conditional_opportunist"),
        "uncertain": sum(1 for r in recs if r.get("sequence_status") in ("candidate_signal", "ambiguous_signal")
                         or (r.get("display_qualifier") is not None and r.get("sequence_status") in _SUPPORTED)),
    })
    pathogens["counts"] = counts


def expected_breadth(depth: float) -> float:
    """Breadth reads at random would give at this mean depth."""
    return 1.0 - math.exp(-depth)
