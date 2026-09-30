"""One sample, end to end: lanes, confirmation, quantification, strains,
evidence, score - into the ``openbiota.mycobiome.v3`` record (spec §13).

Every stage records its own status. A stage that could not run is
``not_assessed`` with a reason; a stage that ran and found nothing is a
result. The two are never the same value.
"""

from __future__ import annotations

import re
import time
from collections.abc import Mapping
from dataclasses import asdict
from pathlib import Path
from typing import Any, Final

from openbiota.mycobiome import (
    eukdetect,
    genome_lane,
    quantify,
    references,
    registry,
    score,
    strains,
)

SCHEMA: Final = "openbiota.mycobiome.v3"
#: Marker-lane support model (EukDetect2's own): accepted as supported when
#: at least this many distinct markers and reads; otherwise provisional.
MARKER_SUPPORT_MARKERS: Final = 3
MARKER_SUPPORT_READS: Final = 8


def _gates(results: dict[str, Any]) -> dict[str, Any]:
    return ((results.get("sequencing_quality") or {}).get("gates") or {})


def _gate_value(results: dict[str, Any], name: str) -> Any:
    for g in _gates(results).get("gates") or []:
        if g.get("name") == name:
            return g.get("value")
    return None


def _organism_names(o: Mapping[str, Any]) -> list[str]:
    """Every name an inventory record answers to: its key, its GTDB name,
    its former name and the names it was listed under before merging."""
    names = [str(o.get("species", "")), str(o.get("gtdb") or ""), str(o.get("formerly") or ""),
             str((o.get("native_ids") or {}).get("scoring") or ""),
             *[str(a) for a in (o.get("aliases") or ())]]
    return [n.replace("_", " ") for n in names if n]


def _bacteria_present(results: dict[str, Any]) -> list[str]:
    inv = results.get("organism_inventory") or {}
    return [n for o in inv.get("organisms") or [] for n in _organism_names(o)]


def pathogen_fungal_signals(results: dict[str, Any]) -> list[dict[str, Any]]:
    """Every fungal target of the pathogen screen that carried at least one
    fragment, with the screen's own status and words. The two modules must
    tell one story: what the screen calls a candidate signal, this page
    shows as the same trace, never as a detection."""
    out: list[dict[str, Any]] = []
    for r in (results.get("pathogens") or {}).get("results") or []:
        tid = str(r.get("target_id") or "")
        if not tid.startswith("fungi."):
            continue
        frags = int(r.get("unique_supporting_fragments") or 0)
        if frags <= 0:
            continue
        name = str(r.get("display_name") or tid.split(".", 1)[1].replace("_", " ").capitalize())
        out.append({
            "target_id": tid, "name": name, "display_status": r.get("display_status"),
            "fragments": frags, "regions": r.get("informative_regions_supported"),
            "identity": r.get("median_alignment_identity"), "statement": r.get("plain_statement") or "",
            "counts_as_pathogen": bool(r.get("counts_as_pathogen")),
        })
    out.sort(key=lambda x: -x["fragments"])
    return out


def _pathogen_fungal_names(results: dict[str, Any]) -> set[str]:
    return {s["name"] for s in pathogen_fungal_signals(results)}


def analyze(*, sample: str, r1: Path, r2: Path | None, results: dict[str, Any], out_dir: Path,
            refs_dir: Path, threads: int = 16, read_sketch: Path | None = None, log=print) -> dict[str, Any]:
    t0 = time.time()
    base = refs_dir / "mycobiome"
    work = out_dir / "mycobiome"
    work.mkdir(parents=True, exist_ok=True)
    lock = references.load_lock(base)
    reg = registry.load()
    eligible = _gates(results).get("usable_nonhost_pairs")
    read_len = _gate_value(results, "read_length") or 150
    rec: dict[str, Any] = {
        "schema_version": SCHEMA, "sample_id": sample, "assay": "stool_shotgun_dna",
        "reference_lock_id": lock.lock_id if lock else None, "input_manifest_id": f"{r1.name}|{r2.name if r2 else ''}",
        "analysis_status": "not_assessed", "stages": {}, "limits": [], "next_steps": [], "exposure_context": [],
        "reference_comparisons": [], "validation_release_id": None,
    }
    if lock is None:
        rec["analysis_status"] = "not_assessed"
        rec["limits"].append("Fungal reference lock not built: run `openbiota mycobiome references prepare`.")
        rec["health_score"] = score.compute([], reason_if_empty="incomplete_analysis", incomplete=True).to_json()
        return rec

    # ---- marker lane -------------------------------------------------------- #
    marker = eukdetect.run_eukdetect(sample=sample, r1=r1, r2=r2, db_dir=base / "eukdetect2", work_dir=work, threads=threads)
    rec["stages"]["marker_lane"] = marker.status if marker else "not_assessed"
    rec["marker_lane"] = marker.to_json() if marker else {"status": "not_assessed", "reason": "EukDetect2 or its database not installed"}
    name_to_tid = {a.species: a.species_taxid for a in lock.assemblies if a.species_taxid}
    extra_species: set[int] = set()
    if marker:
        for h in marker.fungal:
            tid = name_to_tid.get(h.name) or (h.taxid if h.rank == "species" else None)
            if tid:
                extra_species.add(tid)
    for name in _pathogen_fungal_names(results):
        e = reg.lookup(name)
        tid = name_to_tid.get(e.name if e else name)
        if tid:
            extra_species.add(tid)

    # ---- whole-genome lane ---------------------------------------------------- #
    if not genome_lane.available(lock, base):
        rec["stages"]["genome_lane"] = "not_assessed"
        rec["genome_lane"] = {"status": "not_assessed", "reason": "sketch database or aligner not available"}
        genome = None
    else:
        genome = genome_lane.run_genome_lane(
            sample=sample, r1=r1, r2=r2, lock=lock, base=base, work_dir=work, eligible_fragments=int(eligible or 0),
            extra_species_taxids=extra_species, extra_accessions=set(), read_sketch=read_sketch, threads=threads, log=log,
        )
        rec["stages"]["genome_lane"] = genome.status
        rec["genome_lane"] = genome.to_json()
    return _finish(rec, sample=sample, lock=lock, reg=reg, marker=marker, genome=genome, results=results,
                   eligible=eligible, read_len=int(read_len), base=base, work=work, t0=t0, log=log)


def _taxa_records(lock, reg, marker, genome) -> list[dict[str, Any]]:
    """One record per fungal finding, from both lanes, keyed by species."""
    by_acc = lock.by_accession()
    recs: dict[str, dict[str, Any]] = {}
    if genome is not None:
        for t in genome.taxa:
            if t.detection_state == "not_supported" and t.fragments < 2:
                continue
            e = reg.lookup(t.name)
            recs[t.key] = {
                "finding_id": f"{t.key}", "name": t.name, "accepted_name": e.name if e else t.name,
                "rank": t.rank, "species_taxid": t.species_taxid, "genus": t.genus,
                "detection_state": t.detection_state, "support_basis": "whole_genome_competitive",
                "reasons": t.reasons, "genome_lane": t.to_json(), "marker_lane": None,
                "category": e.category if e else "unknown", "category_words": e.category_words if e else "no health evidence",
                "role": e.role if e else "", "origin_interpretation": [e.origin_default] if e else ["unknown"],
                "evidence": [asdict(c) | {"label_words": c.label_words} for c in (e.evidence if e else ())],
                "must_not_infer": e.must_not_infer if e else "", "complex": e.complex if e else None,
                "accessions": t.accessions, "backbone": t.backbone,
                "backbone_organism": by_acc[t.backbone].organism if t.backbone in by_acc else None,
            }
    if marker is not None:
        for h in marker.fungal:
            key = next((k for k, r in recs.items() if r["species_taxid"] == h.taxid or r["name"] == h.name), None)
            if key is not None:
                recs[key]["marker_lane"] = h.to_json()
                continue
            e = reg.lookup(h.name)
            strong = (h.observed_markers or 0) >= MARKER_SUPPORT_MARKERS and h.total_reads >= MARKER_SUPPORT_READS
            recs[f"marker:{h.taxid}"] = {
                "finding_id": f"marker:{h.taxid}", "name": h.name, "accepted_name": e.name if e else h.name,
                "rank": h.rank if h.rank in ("species", "genus") else "higher", "species_taxid": h.taxid if h.rank == "species" else None,
                "genus": h.lineage.split(";")[-2] if ";" in h.lineage else "",
                "detection_state": "supported" if strong and h.rank == "species" else "provisional",
                "support_basis": "marker_lane_only", "reasons": [] if strong else ["below_marker_lane_support_model"],
                "genome_lane": None, "marker_lane": h.to_json(),
                "category": e.category if e else "unknown", "category_words": e.category_words if e else "no health evidence",
                "role": e.role if e else "", "origin_interpretation": [e.origin_default] if e else ["unknown"],
                "evidence": [asdict(c) | {"label_words": c.label_words} for c in (e.evidence if e else ())],
                "must_not_infer": e.must_not_infer if e else "", "complex": e.complex if e else None,
                "accessions": [], "backbone": None, "backbone_organism": None,
            }
    order = {"supported": 0, "ambiguous_complex": 1, "provisional": 2, "trace": 3, "not_supported": 4}
    return sorted(recs.values(), key=lambda r: (order.get(r["detection_state"], 9),
                                                -((r.get("genome_lane") or {}).get("fragments") or 0)))


def _cross_reference_pathogen_screen(taxa: list[dict[str, Any]], reg, results: dict[str, Any]) -> None:
    """Attach the pathogen screen's fungal signals to the matching taxa and
    list the rest as traces under the screen's own status words."""
    signals = pathogen_fungal_signals(results)
    for sig in signals:
        # "Candida krusei (Pichia kudriavzevii)": both names are the same organism.
        names = [n.strip() for n in re.split(r"[()]", sig["name"]) if n.strip()]
        e = next((reg.lookup(n) for n in names if reg.lookup(n)), None)
        accepted = e.name if e else names[0]
        candidates = set(names) | {accepted}
        for r in taxa:
            if r["accepted_name"] in candidates or r["name"] in candidates:
                r["pathogen_screen"] = sig
                break
        else:
            taxa.append({
                "finding_id": f"pathogen:{sig['target_id']}", "name": sig["name"], "accepted_name": accepted,
                "rank": "species", "species_taxid": e.taxid if e else None, "genus": accepted.split(" ")[0],
                "detection_state": "trace", "support_basis": "pathogen_screen_only",
                "reasons": [f"pathogen screen: {sig['display_status']} ({sig['fragments']} fragment{'s' if sig['fragments'] != 1 else ''})"],
                "genome_lane": None, "marker_lane": None, "pathogen_screen": sig,
                "category": e.category if e else "unknown", "category_words": e.category_words if e else "no health evidence",
                "role": e.role if e else "", "origin_interpretation": [e.origin_default] if e else ["unknown"],
                "evidence": [asdict(c) | {"label_words": c.label_words} for c in (e.evidence if e else ())],
                "must_not_infer": e.must_not_infer if e else "", "complex": e.complex if e else None,
                "accessions": [], "backbone": None, "backbone_organism": None,
            })


def _finish(rec, *, sample, lock, reg, marker, genome, results, eligible, read_len, base, work, t0, log):
    log(f"  mycobiome: assembling record for {sample}")
    ledger = genome.ledger if genome is not None else None
    denominator_valid = bool(eligible) and (results.get("sequencing_quality") is not None)
    meas = quantify.measurement(ledger, eligible=eligible, denominator_valid=denominator_valid,
                                lane_status=genome.status if genome else "not_assessed")
    taxa_support = genome.taxa if genome is not None else []
    comp = quantify.composition(taxa_support, fragment_bases=2.0 * read_len, ledger=ledger)
    rec["measurement"] = meas.to_json()
    rec["within_fungi"] = comp.to_json()
    rec["diversity"] = quantify.diversity(taxa_support).to_json()
    taxa = _taxa_records(lock, reg, marker, genome)
    _cross_reference_pathogen_screen(taxa, reg, results)
    rec["pathogen_screen_fungal_signals"] = pathogen_fungal_signals(results)
    rec["searched"] = _searched(lock)

    # ---- mandatory strain attempt for every supported taxon ------------------- #
    strain_recs: list[strains.StrainRecord] = []
    failures: list[str] = []
    eligible_bases = float(eligible or 0) * 2.0 * read_len
    if genome is not None and genome.bam:
        for t in taxa_support:
            if t.detection_state != "supported" or t.rank != "species":
                continue
            frac = (t.fragments / eligible) if eligible else None
            sr = strains.assess(taxon=t, bam=Path(genome.bam), lock=lock, base=base, work_dir=work / "strains",
                                eligible_bases=eligible_bases, fungal_fraction=frac, log=log)
            strain_recs.append(sr)
            if sr.analysis_status == "failed":
                failures.append(f"{t.name}: {sr.error}")
    by_key = {s.taxon_key: s for s in strain_recs}
    for r in taxa:
        s = by_key.get(r["finding_id"])
        r["strain"] = s.to_json() if s else strains.StrainRecord(taxon_key=r["finding_id"], species=r["name"]).to_json()
    required = sum(1 for t in taxa_support if t.detection_state == "supported" and t.rank == "species")
    completed = sum(1 for s in strain_recs if s.analysis_status == "complete")
    resolved = sum(1 for s in strain_recs if s.resolution in ("reference_genotype", "lineage", "reference_equivalence_group", "sample_genotype"))
    rec["strain_analysis"] = {
        "required": True,
        "status": "complete" if required == completed and not failures else ("partial" if strain_recs else ("not_assessed" if genome is None else "complete")),
        "panel_lock_id": lock.lock_id, "taxa_eligible": required, "taxa_assessed": completed,
        "taxa_with_resolved_reference_or_lineage": resolved, "execution_failures": failures, "pairwise_comparisons": [],
        "panels_available": sorted(p.stem.split(".")[0] for p in strains.panels_dir(base).glob("*.panel.json")),
    }
    rec["taxa"] = taxa

    # ---- evidence and score ------------------------------------------------- #
    supported = [
        registry.SupportedTaxon(
            finding_id=r["finding_id"], name=r["accepted_name"], rank=r["rank"], detection_state=r["detection_state"],
            strain_resolution=r["strain"]["resolution"], reference_accessions=tuple(r["strain"]["reference_accessions"]),
            reference_genotypes=tuple(r["strain"]["reference_genotypes"]) if r["strain"]["resolution"] == "reference_genotype" else (),
            compatible_genotypes=tuple(r["strain"]["compatible_genotypes"]),
            compatible_accessions=tuple(a for a in r["strain"]["reference_accessions"]) if r["strain"]["compatible_genotypes"] else (),
            strain_analysis_status=r["strain"]["analysis_status"],
        )
        for r in taxa
    ]
    contribs = registry.evaluate(reg, supported, bacteria_present=_bacteria_present(results))
    {c.rule_id for c in contribs if c.activation_status == "active"}
    scored_fids = {f for c in contribs if c.activation_status == "active" for f in c.finding_ids}
    f_total = ledger.fungal_confident if ledger else 0
    scored_frag = sum((r.get("genome_lane") or {}).get("fragments") or 0 for r in taxa if r["finding_id"] in scored_fids)
    unscored = [r["accepted_name"] for r in taxa if r["detection_state"] in ("supported", "provisional", "ambiguous_complex")
                and r["finding_id"] not in scored_fids]
    unresolved_concerns = [c.rule_id for c in contribs if c.activation_status == "unresolved" and c.direction == "C"]
    alerts = [r["accepted_name"] for r in taxa if r["category"] == "pathogen_alert_route" and r["detection_state"] == "supported"]
    lanes_ok = genome is not None and genome.status == "resolved"
    conf = "insufficient" if not lanes_ok else ("supported" if any(t.detection_state == "supported" for t in taxa_support) else "limited")
    incomplete = (genome is None or genome.status != "resolved" or bool(failures))
    reason = "failed_qc" if _gates(results).get("overall") == "fail" else (
        "incomplete_analysis" if incomplete else ("insufficient_fungal_information" if not taxa else "no_directional_evidence"))
    sc = score.compute(
        contribs, reason_if_empty=reason, analytical_confidence=conf, incomplete=incomplete, unscored_findings=unscored,
        unresolved_concerns=unresolved_concerns, active_alert_ids=alerts, context_missing=[],
        scored_fungal_fragment_fraction=(scored_frag / f_total) if f_total else None,
        strain_assessment_completion={"required": required, "completed": completed, "resolved": resolved,
                                      "completion_fraction": (completed / required) if required else None,
                                      "resolved_fraction": (resolved / required) if required else None},
        policy_lock_id=reg.policy_lock_id,
    )
    rec["health_score"] = sc.to_json()
    rec["context"] = _context(taxa, meas)
    ctx = rec["context"]
    opp_supported = sum(1 for r in taxa if r["category"] in ("commensal_opportunist", "opportunist")
                        and r["detection_state"] == "supported" and r["rank"] == "species")
    sacch = any(r["accepted_name"] == "Saccharomyces cerevisiae" and r["detection_state"] == "supported" for r in taxa)
    rec["myco_score"] = score.myco_score(
        opportunist_share=ctx.get("opportunist_share_of_fungal_fragments"),
        opportunist_fraction_all=ctx.get("opportunist_fraction_of_all_fragments"), opportunist_supported=opp_supported,
        alert_supported=len(alerts), saccharomyces_supported=sacch, benefit=sc.benefit_contribution or 0.0,
        concern=sc.concern_contribution or 0.0, complete=(lanes_ok and not failures),
    )
    rec["colonisation_resistance"] = _colonisation_resistance(results, ctx)
    rec["beneficial_evidence_findings"] = [c.to_json() for c in contribs if c.direction == "B" and c.activation_status != "inactive"]
    rec["potential_concern_findings"] = [c.to_json() for c in contribs if c.direction == "C" and c.activation_status != "inactive"]
    rec["rules_evaluated"] = [c.to_json() for c in contribs]
    rec["analysis_status"] = "complete" if lanes_ok and not failures else ("partial" if lanes_ok else "failed" if genome and genome.status == "failed" else "not_assessed")
    rec["next_steps"] = _next_steps(taxa, sc, alerts)
    rec["limits"] = _limits(marker, genome, meas, lock)
    rec["reference_counts"] = lock.counts
    rec["timings_s"] = {"total": round(time.time() - t0, 1), **(genome.timings_s if genome else {})}
    return rec


#: Groups a reader asks about, and the genera that make them up.
SEARCH_GROUPS: Final = (
    ("Candida and relatives", ("Candida", "Nakaseomyces", "Candidozyma", "Clavispora", "Meyerozyma", "Pichia", "Debaryomyces",
                               "Kluyveromyces", "Cyberlindnera", "Wickerhamomyces", "Yarrowia", "Lodderomyces")),
    ("Saccharomyces and brewing/baking yeasts", ("Saccharomyces", "Brettanomyces", "Torulaspora", "Lachancea", "Kazachstania",
                                                  "Zygosaccharomyces", "Hanseniaspora")),
    ("Malassezia and other skin yeasts", ("Malassezia", "Trichosporon", "Cutaneotrichosporon", "Rhodotorula", "Cryptococcus")),
    ("Aspergillus", ("Aspergillus",)),
    ("Penicillium", ("Penicillium", "Talaromyces")),
    ("Mucor, Rhizopus and relatives", ("Mucor", "Rhizopus", "Rhizomucor", "Lichtheimia", "Cunninghamella")),
    ("Edible mushrooms", ("Agaricus", "Lentinula", "Pleurotus", "Flammulina", "Hericium", "Ganoderma")),
    ("Other moulds", ("Cladosporium", "Alternaria", "Fusarium", "Trichoderma", "Aureobasidium", "Botrytis")),
)


def _searched(lock) -> dict[str, Any]:
    """How many species were searched for, by the groups a reader asks about."""
    by_species = lock.species_assemblies()
    genus_of = {tid: asms[0].genus for tid, asms in by_species.items()}
    groups = []
    covered: set[int] = set()
    for label, genera in SEARCH_GROUPS:
        tids = [tid for tid, g in genus_of.items() if g in genera]
        covered |= set(tids)
        groups.append({"group": label, "species": len(tids), "assemblies": sum(len(by_species[t]) for t in tids)})
    groups.append({"group": "All other fungi", "species": len(by_species) - len(covered),
                   "assemblies": sum(len(a) for t, a in by_species.items() if t not in covered)})
    return {"species_total": len(by_species), "assemblies_total": sum(len(a) for a in by_species.values()),
            "genera_total": len(set(genus_of.values())), "groups": groups,
            "marker_lane": "EukDetect2: single-copy marker genes across its eukaryotic database, filtered to fungi",
            "pathogen_screen": "targeted fungal pathogen catalogue (separate module), cross-referenced here"}


#: Bacteria the review names as suppressing Candida: B. thetaiotaomicron and
#: Blautia producta (LL-37 via HIF-1α, mouse; Fan et al. 2015 Nat Med) and
#: Lactobacillus rhamnosus (metabolites and nutrient competition, mouse and
#: in vitro). Checked against the sample's own bacterial inventory.
COLONISATION_RESISTANCE_BACTERIA: Final = (
    ("Bacteroides thetaiotaomicron", "LL-37 induction via HIF-1\u03b1; resists C. albicans colonisation (mouse)"),
    ("Blautia producta", "LL-37 induction via HIF-1\u03b1; resists C. albicans colonisation (mouse)"),
    ("Lacticaseibacillus rhamnosus", "reduces C. albicans pathogenicity by metabolites and nutrient competition (mouse, in vitro)"),
)
COLONISATION_RESISTANCE_SOURCE: Final = "Huang et al. 2024, Gut Microbes 16:2440111, section 4.3 (refs 148, 149)"


def _level_percentile_of(rec: Mapping[str, Any]) -> float | None:
    from openbiota.inventory import level_percentile_of

    return level_percentile_of(rec)


def _colonisation_resistance(results: dict[str, Any], ctx: dict[str, Any]) -> dict[str, Any]:
    """Which Candida-suppressing bacteria this sample carries, from the
    bacterial inventory. Context for an opportunist finding; never a score input."""
    inv: dict[str, Any] = {}
    for o in (results.get("organism_inventory") or {}).get("organisms") or []:
        for n in _organism_names(o):
            inv.setdefault(n, o)
    alias = {"Lacticaseibacillus rhamnosus": ("Lactobacillus rhamnosus",)}
    rows = []
    for name, role in COLONISATION_RESISTANCE_BACTERIA:
        hit = inv.get(name) or next((inv[a] for a in alias.get(name, ()) if a in inv), None)
        rows.append({"bacterium": name, "role": role, "present": hit is not None,
                     "percent": (hit.get("percent") if hit else None),
                     "percentile": (_level_percentile_of(hit) if hit else None)})
    return {"source": COLONISATION_RESISTANCE_SOURCE, "bacteria": rows,
            "relevant": bool(ctx.get("opportunist_taxa"))}


#: Published context for a reader's "is this normal?" - the HMP mycobiome
#: survey (Nash et al. 2017, Microbiome 5:153): 317 stool samples from 147
#: healthy US adults. Prevalence is ITS2 amplicon (a different assay from
#: this one); the shotgun figure is method-specific. Context, never a target.
HMP_SOURCE: Final = "https://link.springer.com/article/10.1186/s40168-017-0373-4"
HMP_SHOTGUN_FUNGAL_FRACTION: Final = 0.0001      # ~0.01% of >27 billion shotgun reads mapped to fungi
HMP_GENUS_PREVALENCE: Final = (                 # genus, % of samples (ITS2), what it usually means
    ("Saccharomyces", 96.8, "bread, beer, wine, probiotics"),
    ("Malassezia", 88.3, "skin yeast; stool carriage usually from skin"),
    ("Candida", 80.8, "common resident; C. albicans most often"),
)
#: Genera reported as common in healthy stool without a single prevalence
#: figure to quote (Nash 2017; Raimondi 2019; Hallen-Adams & Suhr 2017; Huang
#: et al. 2024). Each says what carrying it usually means: a column reading
#: "commonly reported" twice over tells a reader nothing they did not already
#: have from the row's presence in the table.
HMP_GENUS_COMMON: Final = (
    ("Penicillium", "mould on cheese and cured food; also inhaled from air"),
    ("Debaryomyces", "salt-tolerant yeast of cheese rind and cured meat"),
    ("Cladosporium", "one of the commonest airborne moulds; swallowed daily"),
    ("Aspergillus", "airborne mould; stool reads usually mean exposure, not colonisation"),
    ("Pichia", "fermentation yeast; competes with Candida in the gut"),
    ("Kluyveromyces", "dairy yeast of kefir and some cheeses"),
    ("Cyberlindnera", "fermentation yeast, often from plant food"),
    ("Geotrichum", "rind yeast of soft cheese; a normal transient"),
    ("Rhodotorula", "pink environmental yeast of damp surfaces and produce"),
)

CLASS_GROUPS: Final = (
    ("food_and_drink", "food & drink yeasts", ("food_associated", "studied_probiotic_species")),
    ("opportunist", "opportunistic yeasts", ("commensal_opportunist", "opportunist")),
    ("beneficial_preclinical", "preclinical beneficial", ("beneficial_preclinical",)),
    ("environmental", "environmental & food moulds", ("environmental_or_food",)),
    ("pathogen_route", "pathogen alert route", ("pathogen_alert_route",)),
    ("unknown", "no health evidence", ("unknown",)),
)


def _context(taxa: list[dict[str, Any]], meas) -> dict[str, Any]:
    """Healthy-adult context and the class make-up of the fungal DNA found."""
    listed = [t for t in taxa if t["detection_state"] in ("supported", "ambiguous_complex", "provisional", "trace")]
    # A former Candida species is written "[Candida] boidinii". Matching the
    # bracketed form literally left the Candida row reading "not detected" on
    # a page that listed a [Candida] trace three inches above it.
    def _genus_key(t: dict[str, Any]) -> str:
        return str(t.get("genus") or t["accepted_name"].split(" ")[0]).strip("[]")

    genera_found = {_genus_key(t): t["detection_state"] for t in listed}
    rows = []
    for genus, prev, meaning in HMP_GENUS_PREVALENCE:
        rows.append({"genus": genus, "hmp_prevalence_percent": prev, "meaning": meaning,
                     "found_here": genera_found.get(genus, "not_detected")})
    for genus, meaning in HMP_GENUS_COMMON:
        rows.append({"genus": genus, "hmp_prevalence_percent": None, "meaning": meaning,
                     "found_here": genera_found.get(genus, "not_detected")})
    frac = meas.fungal_fragment_fraction
    if frac is None:
        position = "not measured"
    elif frac < HMP_SHOTGUN_FUNGAL_FRACTION / 10:
        position = "below the HMP figure by more than tenfold"
    elif frac > HMP_SHOTGUN_FUNGAL_FRACTION * 10:
        position = "above the HMP figure by more than tenfold"
    else:
        position = "within tenfold of the HMP figure"
    classes = []
    counted = [t for t in listed if t["detection_state"] != "trace"]

    def _frags(t: dict[str, Any]) -> int:
        return (t.get("genome_lane") or {}).get("fragments") or (t.get("pathogen_screen") or {}).get("fragments") or 0

    counted_total = sum(_frags(t) for t in counted)
    for key, label, cats in CLASS_GROUPS:
        members = [t for t in counted if t["category"] in cats]
        frags = sum(_frags(t) for t in members)
        classes.append({"key": key, "label": label, "taxa": [t["accepted_name"] for t in members], "fragments": frags,
                        "share_of_fungal_fragments": (frags / counted_total) if counted_total else None})
    opp = next(c for c in classes if c["key"] == "opportunist")
    alert = next(c for c in classes if c["key"] == "pathogen_route")
    eligible = meas.eligible_fragments or 0
    opp_all = (opp["fragments"] / eligible) if eligible else None
    return {
        "opportunist_fraction_of_all_fragments": opp_all,
        "opportunist_fragments": opp["fragments"],
        "counted_fungal_fragments": counted_total,
        "source": HMP_SOURCE,
        "hmp_shotgun_fungal_fraction": HMP_SHOTGUN_FUNGAL_FRACTION,
        "sample_fungal_fraction": frac,
        "position_versus_hmp": position,
        "genus_prevalence": rows,
        "note": ("HMP prevalence is ITS2 amplicon on 317 stool samples from 147 healthy US adults; the shotgun read "
                 "fraction is method-specific. Neither is a healthy range or a target: less fungal DNA is not better, "
                 "and more is not worse."),
        "class_composition": classes,
        "opportunist_share_of_fungal_fragments": opp["share_of_fungal_fragments"],
        "opportunist_taxa": opp["taxa"],
        "alert_route_taxa": alert["taxa"],
    }


def _next_steps(taxa, sc, alerts) -> list[dict[str, str]]:
    """Result-specific follow-up (spec §11), attached to findings, not the gauge."""
    out: list[dict[str, str]] = []
    sup = [r for r in taxa if r["detection_state"] == "supported"]
    if not sup:
        out.append({"trigger": "no_supported_calls", "text": (
            "No fungus reached the support rule. The whole-genome search covered every admitted fungal assembly "
            "and the marker search every eukaryotic marker in its database; at this depth a fungus below roughly "
            "0.001% of fragments can be missed. If fungal symptoms persist, a targeted fungal assay or a repeat "
            "sample answers what low-depth shotgun cannot.")})
    if any(r["category"] in ("studied_probiotic_species", "food_associated") for r in sup):
        out.append({"trigger": "food_or_probiotic_signal", "text": (
            "Food- or probiotic-associated yeasts were found. Recording yeast probiotic and fermented-food intake and "
            "the time since the last dose separates a transient passenger from a persistent resident; a repeat "
            "sample after a break is the standard way to tell. Do not stop a beneficial therapy on a sequence result alone.")})
    if any(r["category"] == "beneficial_preclinical" for r in sup):
        out.append({"trigger": "preclinical_beneficial_species", "text": (
            "A species with preclinical beneficial evidence was found, but the studied strain was not established in "
            "this sample. The research is shown with that limit; it is not a supplementation recommendation.")})
    if sc.unresolved_concerns:
        out.append({"trigger": "strain_unresolved_concern", "text": (
            "A concern that depends on the exact strain could not be resolved at this depth. The organism's page "
            "lists the loci searched; isolate sequencing or a deeper run would settle it. Species presence alone is "
            "not the studied phenotype.")})
    if alerts:
        out.append({"trigger": "pathogen_alert_route", "text": (
            f"{', '.join(alerts)} is routed to the pathogen alert policy. This stool module does not establish invasive "
            "or systemic fungal infection; the appropriate clinical evaluation does.")})
    out.append({"trigger": "always", "text": (
        "This module reads stool DNA. It does not establish small-intestinal fungal overgrowth, mucosal invasion, "
        "fungemia, viability, growth phase or toxin production. No eradication target for all fungi follows from "
        "any result here, and no antifungal is selected by relative abundance.")})
    return out


def _limits(marker, genome, meas, lock) -> list[str]:
    lim = [
        "Fungal DNA is not living biomass: genome size, cell-wall lysis, ploidy and dietary DNA all shape the signal.",
        "Within-fungi percentages and the total-sample fungal signal use different denominators and are shown apart.",
        f"Reference: {lock.counts.get('assemblies_admitted_fungal', 0):,} fungal assemblies covering "
        f"{lock.counts.get('species_admitted', 0):,} species; a fungus outside them is found only by the marker lane, "
        "and one outside both catalogues is not found at all.",
        "The support rule (20 fragments over 3 separated regions) is a locked starting value from the specification's "
        "calibration grid; qualification on mock communities is pending and is recorded as such.",
        "Stool does not fully represent the intestinal mucosa; a mucosal fungus can be absent from stool.",
    ]
    if marker is None:
        lim.append("The marker lane (EukDetect2) was not run.")
    if genome is None or genome.status != "resolved":
        lim.append("The whole-genome lane did not complete; detections here are not exhaustive.")
    if meas.quantification_status == "detected_below_quantification_limit":
        lim.append("Fungal fragments are below the quantification limit; the count is exact, the percentage imprecise.")
    return lim


__all__ = ["SCHEMA", "analyze"]
