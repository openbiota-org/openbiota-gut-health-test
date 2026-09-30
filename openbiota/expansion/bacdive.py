"""BacDive enrichment: culture-collection identities and reference-strain traits.

Spec 0.8.4 §6. BacDive (DSMZ) is strain and phenotype metadata, not a
classifier: what it adds to a detected species is its culture-collection
identities, linked genome accessions, oxygen requirement, temperature
range, substrates used and products formed, and where type strains were
isolated. Every trait is a trait of a *reference strain* and is labelled
so; it is never asserted of the patient's organism.

API v2, freely accessible without registration since February 2026.
`/v2/taxon/{genus}/{species}` returns BacDive ids (paginated);
`/v2/fetch/{id1;id2;...}` returns records, at most 100 ids per call. Raw
responses are cached under refs/bacdive/cache with timestamps; requests
retry with backoff; a species with no record is recorded as
`no_matching_record`.

Licence: CC BY 4.0; DSMZ asks to be contacted for commercial use. The
terms travel with the block written to results.json.
"""

from __future__ import annotations

import contextlib
import datetime as dt
import hashlib
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any, Final

BASE: Final = "https://api.bacdive.dsmz.de/v2"
CACHE: Final = Path("refs/bacdive/cache")
BATCH: Final = 100
MAX_STRAINS_PER_SPECIES: Final = 12
LICENCE: Final = ("BacDive data are CC BY 4.0 (DSMZ). DSMZ asks to be contacted for commercial use. "
                  "Each record carries its own DOI for citation.")
USER_AGENT: Final = "openbiota/0.8.4 (research; https://openbiota.org)"


def _get(url: str, *, retries: int = 4) -> Any:
    key = hashlib.sha256(url.encode()).hexdigest()[:24]
    CACHE.mkdir(parents=True, exist_ok=True)
    cached = CACHE / f"{key}.json"
    if cached.is_file():
        try:
            return json.loads(cached.read_text())["body"]
        except (json.JSONDecodeError, KeyError):
            pass
    delay = 1.5
    last: Exception | None = None
    for _ in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT, "Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=60) as resp:
                body = json.loads(resp.read().decode("utf-8"))
            cached.write_text(json.dumps({"url": url, "fetched_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
                                          "body": body}))
            return body
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return None
            last = exc
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError) as exc:
            last = exc
        time.sleep(delay)
        delay *= 2
    raise RuntimeError(f"BacDive request failed after {retries} attempts: {url}: {last}")


def taxon_ids(genus: str, species: str) -> list[int]:
    """Every BacDive id for a species, following pagination."""
    url = f"{BASE}/taxon/{urllib.parse.quote(genus)}/{urllib.parse.quote(species)}"
    ids: list[int] = []
    while url:
        body = _get(url)
        if not body:
            break
        ids.extend(int(x) for x in body.get("results") or [])
        url = body.get("next")
    return ids


def fetch(ids: list[int]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for i in range(0, len(ids), BATCH):
        chunk = ids[i:i + BATCH]
        body = _get(f"{BASE}/fetch/{';'.join(str(x) for x in chunk)}")
        if not body:
            continue
        results = body.get("results")
        if isinstance(results, dict):
            out.extend(results.values())
        elif isinstance(results, list):
            out.extend(results)
    return out


def _first(x: Any) -> Any:
    if isinstance(x, list):
        return x[0] if x else None
    return x


def summarise(records: list[dict[str, Any]]) -> dict[str, Any]:
    """The traits a reader can use, aggregated over a species' reference strains."""
    culture_ids: list[str] = []
    genome_accessions: list[str] = []
    oxygen: dict[str, int] = {}
    temps: list[float] = []
    substrates_used: dict[str, int] = {}
    products: dict[str, int] = {}
    isolation: dict[str, int] = {}
    dois: list[str] = []
    type_strains = 0
    for rec in records:
        gen = rec.get("General") or {}
        if gen.get("DSM-Number"):
            culture_ids.append(f"DSM {gen['DSM-Number']}")
        if gen.get("doi"):
            dois.append(str(gen["doi"]))
        tax = rec.get("Name and taxonomic classification") or {}
        if str(tax.get("type strain", "")).lower() in ("yes", "true"):
            type_strains += 1
        seq = rec.get("Sequence information") or {}
        for g in (seq.get("Genome sequences") if isinstance(seq.get("Genome sequences"), list) else [seq.get("Genome sequences")]):
            if isinstance(g, dict) and g.get("INSDC accession"):
                genome_accessions.append(str(g["INSDC accession"]))
        phys = rec.get("Physiology and metabolism") or {}
        for o in (phys.get("oxygen tolerance") if isinstance(phys.get("oxygen tolerance"), list) else [phys.get("oxygen tolerance")]):
            if isinstance(o, dict) and o.get("oxygen tolerance"):
                k = str(o["oxygen tolerance"])
                oxygen[k] = oxygen.get(k, 0) + 1
        for m in (phys.get("metabolite utilization") or []):
            if isinstance(m, dict) and m.get("utilization activity") == "+" and m.get("metabolite"):
                k = str(m["metabolite"])
                substrates_used[k] = substrates_used.get(k, 0) + 1
        for m in (phys.get("metabolite production") or []):
            if isinstance(m, dict) and str(m.get("production", "")).lower() in ("yes", "+") and m.get("metabolite"):
                k = str(m["metabolite"])
                products[k] = products.get(k, 0) + 1
        cult = rec.get("Culture and growth conditions") or {}
        for t in (cult.get("culture temp") if isinstance(cult.get("culture temp"), list) else [cult.get("culture temp")]):
            if isinstance(t, dict):
                with contextlib.suppress(ValueError):
                    temps.append(float(str(t.get("temperature", "")).split("-")[0]))
        iso = rec.get("Isolation, sampling and environmental information") or {}
        for i in (iso.get("isolation source categories") or []):
            if isinstance(i, dict):
                k = str(i.get("Cat3") or i.get("Cat2") or i.get("Cat1") or "").lstrip("#")
                if k:
                    isolation[k] = isolation.get(k, 0) + 1
    def top(d: dict[str, int], n: int) -> list[str]:
        return [k for k, _ in sorted(d.items(), key=lambda kv: -kv[1])[:n]]
    return {
        "n_strains": len(records), "n_type_strains": type_strains,
        "culture_collection_ids": sorted(set(culture_ids))[:MAX_STRAINS_PER_SPECIES],
        "genome_accessions": sorted(set(genome_accessions))[:MAX_STRAINS_PER_SPECIES],
        "oxygen_requirement": top(oxygen, 2),
        "growth_temperature_c": (round(min(temps), 1), round(max(temps), 1)) if temps else None,
        "substrates_used": top(substrates_used, 12), "products": top(products, 8),
        "isolation_sources": top(isolation, 4), "dois": dois[:MAX_STRAINS_PER_SPECIES],
        "caveat": "traits measured in reference strains; not asserted of the organism in this sample",
    }


def enrich(species_names: list[str], *, budget_s: float = 240.0) -> dict[str, Any]:
    """BacDive summaries for named species, within a time budget.

    Species are taken in the order given (callers pass most abundant
    first). Unnamed clusters have nothing to look up and are skipped.
    """
    t0 = time.monotonic()
    out: dict[str, Any] = {}
    n_hit = n_none = 0
    for name in species_names:
        if time.monotonic() - t0 > budget_s:
            out[name] = {"status": "not_assessed", "reason": "time budget reached"}
            continue
        parts = name.replace("_", " ").split()
        if len(parts) < 2 or not parts[1][:1].islower() or parts[1].startswith("sp"):
            out[name] = {"status": "no_matching_record", "reason": "no binomial to look up"}
            n_none += 1
            continue
        genus = parts[0].split("_")[0]  # GTDB suffixes (_A) are not in LPSN names
        try:
            ids = taxon_ids(genus, parts[1])
            if not ids:
                out[name] = {"status": "no_matching_record"}
                n_none += 1
                continue
            records = fetch(ids[:MAX_STRAINS_PER_SPECIES])
            out[name] = {"status": "found", "bacdive_ids": ids[:MAX_STRAINS_PER_SPECIES], "n_ids_total": len(ids),
                         **summarise(records)}
            n_hit += 1
        except Exception as exc:  # noqa: BLE001 - one failed lookup must not stop the rest
            out[name] = {"status": "failed", "reason": f"{type(exc).__name__}: {str(exc)[:160]}"}
    return {
        "source": "BacDive API v2 (api.bacdive.dsmz.de)", "licence": LICENCE,
        "fetched_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "n_species_queried": len(species_names), "n_found": n_hit, "n_no_record": n_none,
        "elapsed_s": round(time.monotonic() - t0, 1), "species": out,
        "meaning": ("Culture-collection identities, linked genome accessions and reference-strain traits for the named "
                    "species detected. A trait measured in a reference strain is not automatically a trait of the "
                    "strain in this sample."),
    }
