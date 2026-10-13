#!/usr/bin/env python3
"""Turn the research atlas into options the report can put in front of a reader.

    scripts/import_atlas_evidence.py [--zip specs/research/OPENBIOTA_INTERVENTION_DATA.zip]
                                     [--out interventions/atlas_options.yaml] [--report]

The atlas (``specs/research/OPENBIOTA_INTERVENTION_DATA.zip``) holds 87,010
preserved research records. Most of them are database rows, alias tables,
sample indexes and source catalogues: context, not options. About nine
hundred of them are curated observations in exactly the shape an option
needs - an exact material, what it was tested against, what happened, in
what setting, by what route, with its boundary and its citation - and those
are what this reads.

Nothing is filtered on study design. A mouse study and a dish are imported
on the same footing as a randomised trial, because that is where most of
what is known about botanicals, enzymes and antimicrobial food compounds
actually lives. What each row carries instead is a *setting* and a
*delivery route*, and :mod:`openbiota.delivery` decides from those two what
has to be said beside the finding. The only thing that disqualifies a row
is a material nobody swallows.

Targets are stored as the atlas names them. Matching them to the organisms
and panels of a particular report is the runtime's job
(:mod:`openbiota.solutions`), so a new crosswalk does not need a re-import.
"""

from __future__ import annotations

import argparse
import collections
import html
import json
import re
import sqlite3
import sys
import tempfile
import zipfile
from pathlib import Path
from typing import Any

import yaml

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from openbiota import delivery as DEL  # noqa: E402

DEFAULT_ZIP = REPO / "specs" / "research" / "OPENBIOTA_INTERVENTION_DATA.zip"
DEFAULT_OUT = REPO / "interventions" / "atlas_options.yaml"

#: The tables that hold one observation per row: a material, a target and
#: what happened. Everything else in the package is context.
EVIDENCE_TABLES: tuple[str, ...] = (
    "T038", "T039", "T040", "T041", "T042", "T043", "T045", "T046", "T047",
    "T048", "T049", "T050", "T051", "T052", "T053", "T054", "T055", "T056",
    "T057", "T060", "T061", "T062", "T063", "T064",
)
#: Measured ecological chains: material -> first effect -> downstream effect.
#: They answer "and what else would that do", which nothing else here does.
CASCADE_TABLES: tuple[str, ...] = ("T081", "T082", "T083", "T085", "T086")

# --------------------------------------------------------------------------- #
# reading a cell
# --------------------------------------------------------------------------- #

_LINK = re.compile(r"\[([^\]]*)\]\(([^)]*)\)")
_BOLD = re.compile(r"\*\*(.+?)\*\*")
_ITAL = re.compile(r"(?<!\*)\*(?!\*)(.+?)(?<!\*)\*(?!\*)")
_WS = re.compile(r"\s+")


def text_of(cell: Any) -> str:
    """A cell as plain prose: markdown unwrapped, entities decoded."""
    s = html.unescape(str(cell or ""))
    s = s.replace("<br>", " ").replace("<br/>", " ")
    s = _LINK.sub(lambda m: m.group(1) or m.group(2), s)
    s = _BOLD.sub(r"\1", s)
    s = _ITAL.sub(r"\1", s)
    return _WS.sub(" ", s).strip()


def urls_of(cell: Any) -> list[str]:
    raw = html.unescape(str(cell or ""))
    out = [m.group(2) for m in _LINK.finditer(raw) if m.group(2).startswith("http")]
    out += re.findall(r"https?://[^\s;,)\]]+", _LINK.sub("", raw))
    seen: list[str] = []
    for u in out:
        u = u.rstrip(".,;")
        if u not in seen:
            seen.append(u)
    return seen


def pmid_of(blob: str) -> str | None:
    m = re.search(r"pubmed\.ncbi\.nlm\.nih\.gov/(\d+)", blob) or re.search(r"PMID[:\s*]*(\d{6,9})", blob)
    return m.group(1) if m else None


def doi_of(blob: str) -> str | None:
    m = re.search(r"(10\.\d{4,9}/[^\s\"'<>)\]]+)", blob)
    return m.group(1).rstrip(".,;") if m else None


# --------------------------------------------------------------------------- #
# what kind of thing is it, how does it get there, and where was it seen
# --------------------------------------------------------------------------- #

KIND_PATTERNS: tuple[tuple[str, str], ...] = (
    # Every pattern is a prefix match on purpose: the sources write
    # "oligosaccharides", "extracts", "inhibited". A trailing word boundary
    # after a prefix matches none of them, so there is none here, and the
    # alternatives that would be ambiguous as prefixes carry their own.
    ("procedure", r"\b(?:transplant|fmt\b|colonoscop|endoscop|surgery|faecal microbiota transplant)"),
    ("clinical_route", r"\b(?:referral|specialist|gastroenterolog|dermatolog|clinical evaluation|"
                       r"screening programme|guideline route)"),
    ("spore_probiotic", r"\b(?:bacillus |spore|clostridium butyricum|saccharomyces boulardii)"),
    ("probiotic", r"\b(?:lactobacill|bifidobacter|lacticaseibacill|limosilactobacill|lactiplantibacill|"
                  r"streptococcus thermophilus|akkermansia|probiotic|synbiotic|postbiotic|"
                  r"live biotherapeutic|heat-?inactivated|pasteuri[sz]ed|escherichia coli nissle|"
                  r"\becn\b|consortium|defined community)"),
    # Named drugs only. "antibiotic", "antifungal" and "prescription" used to
    # be in here and dragged in anything whose study happened to mention one -
    # which is how a two-strain probiotic came to be labelled a medication.
    ("drug", r"\b(?:rifaximin|vancomycin|metronidazole|fidaxomicin|mesalazine|statin|metformin|"
             r"omeprazole|fluconazole|nystatin|ivermectin|albendazole|neomycin|azithromycin|"
             r"doxycycline|amoxicillin|clarithromycin|levodopa|digoxin|triple therapy|"
             r"proton pump inhibitor)"),
    ("essential_oil", r"\b(?:essential oil|oregano oil|thyme oil|tea tree|carvacrol|thymol|eugenol|"
                      r"peppermint oil|cinnamaldehyde|limonene|oil of )"),
    ("enzyme", r"\b(?:enzyme|lactase|protease|amylase|lysozyme|dnase|serrapeptase|nattokinase|"
               r"alpha-?galactosidase|glucanase|peptidase|lipase)"),
    ("polyphenol", r"\b(?:polyphenol|flavonoid|flavanol|catechin|epicatechin|quercetin|resveratrol|"
                   r"anthocyan|ellagitannin|urolithin|curcumin|berberine|proanthocyanidin|tannin|"
                   r"isoflavone|daidzein|genistein)"),
    ("herb", r"\b(?:extract|botanical|rosemary|oregano|garlic|allicin|ginger|wormwood|artemisia|"
             r"neem|clove|cinnamon|turmeric|green tea|grape seed|pomegranate|olive leaf|barberry|"
             r"goldenseal|liquorice|licorice|aloe|peppermint|fennel|thyme|sage\b|herbal)"),
    ("prebiotic", r"\b(?:inulin|fructo-?oligosacchar|\bfos\b|galacto-?oligosacchar|\bgos\b|"
                  r"xylo-?oligosacchar|\bxos\b|resistant starch|resistant dextrin|arabinoxylan|"
                  r"beta-?glucan|psyllium|guar gum|\bphgg\b|pectin|prebiotic|\bhmo\b|"
                  r"human milk oligosacchar|2'-?fl\b|\blnnt\b|polydextrose|lactulose|oligosacchar|"
                  r"maltodextrin|glucomannan|mucin\b)"),
    ("fibre", r"\b(?:fibre|fiber|cellulose|bran\b|wholegrain|whole grain|legume|oat\b|barley)"),
    ("mineral", r"\b(?:zinc |iron |calcium |magnesium |selenium|bismuth|lactoferrin)"),
    ("food", r"\b(?:yogurt|yoghurt|kefir|kimchi|sauerkraut|miso\b|natto\b|tempeh|cocoa|chocolate|"
             r"walnut|almond|berry|berries|fruit|vegetable|fish oil|olive oil|fermented food|juice|"
             r"milk\b|cheese|bread|coffee|\btea\b|honey|nuts\b|egg\b|meat\b)"),
    ("diet", r"\b(?:diet|fodmap|mediterranean|enteral nutrition|feeding|eating pattern|fasting|\bcded\b)"),
    ("lifestyle", r"\b(?:exercise|sleep|stress |physical activity|smoking)"),
    ("compound", r"\b(?:acid\b|vitamin|n-?acetylcysteine|\bnac\b|butyrate|glutamine|collagen|colostrum|"
                 r"immunoglobulin|melatonin|creatine|supplement|glycerol|reuterin|bacteriocin|"
                 r"microcin|peptide)"),
)

SETTING_PATTERNS: tuple[tuple[str, str], ...] = (
    ("human_trial", r"\b(?:randomi[sz]|\brct\b|placebo|double-?blind|crossover trial|clinical trial|"
                    r"phase [123]|human (?:randomi|trial|challenge|study)|participants|patients were|"
                    r"\bhr\b pilot|adults? (?:were|across)|volunteers)"),
    ("human_observational", r"\b(?:human (?:cohort|observational|fecal|faecal|adults?|volunteers?)|cohort\b|"
                            r"case-?series|case report|observational|cross-?sectional|epidemiolog)"),
    ("animal", r"\b(?:mice\b|mouse\b|murine|\brats?\b|germ-?free|gnotobiotic|piglet|calf\b|calves|"
               r"canine|\bdogs?\b|zebrafish|animal (?:model|study|community)|veterinary)"),
    ("ex_vivo", r"\b(?:ex vivo|faecal slurry|fecal slurry|batch fermentation|simulated (?:colon|gut|digestion)|"
                r"\bshime\b|continuous culture|human tissue|biops|organoid)"),
    ("in_vitro", r"\b(?:in ?vitro|isolate culture|co-?culture|cultur|broth|agar\b|\bmic\b|"
                 r"minimum inhibitory|cell line|caco-?2|ht-?29|biofilm assay|laborator|"
                 r"flow cell|microplate|assay)"),
)

DELIVERY_PATTERNS: tuple[tuple[str, str], ...] = (
    ("not_ingestible", r"\b(?:topical|intravenous|\biv\b|injection|intraperitoneal|inhal|nasal spray|"
                       r"subcutaneous|catheter|wound|skin application)"),
    ("enteric_required", r"\b(?:enteric|delayed-?release|acid-?protect|encapsulat|gastro-?resistant)"),
    ("oral_tested", r"\b(?:oral|by mouth|ingest|\bfed\b|feeding|dietary exposure|capsule|drink|"
                    r"consum|supplementation|administered|\bORAL-[A-Z]\b|tablet|sachet)"),
    ("absorbed_upstream", r"\b(?:absorbed|bioavailab|plasma|serum concentration|systemic)"),
)

#: Direction of the observed effect on the named target.
DIRECTION_PATTERNS: tuple[tuple[str, str], ...] = (
    # Null is tested first on purpose: "no significant reduction" contains
    # "reduc", and reading it as a decrease would invert the finding.
    ("null", r"\b(?:no (?:significant |measurable |detectable |added )?(?:change|improvement|difference|"
             r"effect|benefit|reduction|increase|advantage)|did not (?:significantly )?(?:change|improve|"
             r"differ|reduce|increase|alter|affect|eradicate)|was not (?:significant|superior|improved)|"
             r"null\b|unchanged|failed to|non-?significant|no improvement)"),
    ("decrease", r"\b(?:inhibit|reduc|decreas|lower|suppress|decolon|eradicat|clearance|cleared|kill|"
                 r"cidal\b|antibacterial|antimicrobial|antifungal|deplet|fewer|declin|anti-?adhesion|"
                 r"displac|outcompet|loss of|against\b|lysis|disrupt)"),
    ("increase", r"\b(?:increas|rais|enhanc|promot|stimulat|enrich|restor|higher|greater|growth of|"
                 r"bloom|expand|support(?:ed|s)? (?:the )?growth|bifidogenic|improv|grew\b|positive:)"),
)

ROLE_HINTS: dict[str, tuple[str, ...]] = {
    "material": ("exact ingredient", "exact material", "exact strain", "material", "exact intervention",
                 "exact candidate", "exact studied material", "exact substrate", "exact preparation",
                 "exact target and preparation", "intervention", "downloaded source file",
                 "exact material or formulation", "exact material / intervention", "exact material and biological",
                 "exact material; model/route", "exact material / initial perturbation",
                 "id / exact material", "material / research topic", "material / preparation",
                 "condition 1 / material", "exact ingredient/preparation", "exact ingredient / preparation"),
    "target": ("target", "organism", "report-relevant target", "exact source strain label",
               "target actually measured", "target / endpoint family", "target or endpoint",
               "target and observed role", "target and observed direction", "target / role / observed result",
               "target / setting / observed result", "model and target", "target/study",
               "exact candidate / action target", "exact ingredient / proposed target",
               "taxon and source direction"),
    "effect": ("observed", "actual observed effect", "measured endpoint", "finding", "empirical finding",
               "what was observed", "observed result", "observed endpoint", "observed direction",
               "measured relationship", "evidence and observed direction", "first effect",
               "positive or null result", "verified primary evidence", "edge 1", "conditional chain",
               "intervention → target / observed direction", "input → mediator → output",
               "material → direct target", "potential chain"),
    "setting": ("model", "evidence", "setting", "evidence model", "model and exposure", "model / delivery",
                "evidence/model and delivery", "evidence and context", "evidence/delivery",
                "source host category", "source host / actual method", "chain evidence",
                "end-to-end evidence and endpoint", "e2e measured", "evidence and exposure"),
    "delivery": ("delivery", "delivery evidence", "delivery / preparation", "delivery and commercial identity",
                 "setting and delivery", "delivery and individual-matching interpretation", "live at consumption"),
    "boundary": ("boundary", "limits", "tradeoff", "scope", "interpretation", "breadth", "claim boundary",
                 "limits and interpretation", "limits / collateral", "commercial identity",
                 "boundary or competing finding", "modifiers", "baseline", "existing-intake context",
                 "required interpretation boundary", "permitted action-card claim", "resolution",
                 "edge 3", "edge 2", "downstream"),
    "source": ("primary source", "primary research", "primary publication", "primary doi", "primary study",
               "sources", "source", "primary identifier", "primary or operator source", "archive"),
}


def role_of(header: str) -> str | None:
    h = header.strip().lower()
    if h in {"id", "id / report trigger", "node", "field", "resource"}:
        return "id"
    # Effect is matched before setting on purpose. "Verified primary evidence
    # and endpoint" is the column that states the result, but it contains the
    # word "evidence", so a setting-first order claimed it and left the result
    # blank - which is how rows arrived on the page with nothing in the
    # "what was seen" column.
    for role in ("source", "delivery", "boundary", "effect", "setting", "target", "material"):
        if any(hint in h for hint in ROLE_HINTS[role]):
            return role
    return None


def first_match(patterns: tuple[tuple[str, str], ...], blob: str, default: str) -> str:
    low = blob.lower()
    for name, pattern in patterns:
        if re.search(pattern, low):
            return name
    return default


_TRIGGER_PREFIX = re.compile(r"^\s*(?:report trigger|trigger|finding)\s*:\s*[^;]{0,120};\s*", re.I)
_CATEGORY_SUFFIX = re.compile(r"[;.]\s*categor(?:y|ies)\s*:\s*[^;.]*\.?\s*$", re.I)


def clean_material(cell: str) -> str:
    """The thing itself, without the bookkeeping the atlas wrapped it in.

    Some tables key a row by the report finding it answers, so the material
    cell reads "Report trigger: high BCAA biosynthesis potential; Berberine".
    The reader wants the berberine.
    """
    out = _TRIGGER_PREFIX.sub("", str(cell or "").strip())
    out = _CATEGORY_SUFFIX.sub("", out).strip(" ;.,")
    return out or str(cell or "").strip()


def classify_kind(material: str, whole: str) -> str:
    """What sort of thing it is, judged on its own name first.

    Reading the whole row put an iron-metabolism aside into the kind of a
    probiotic strain; the material's own name decides, and the rest of the
    row is consulted only when the name alone says nothing.
    """
    named = first_match(KIND_PATTERNS, material, "")
    return named or first_match(KIND_PATTERNS, whole, "compound")


def classify_setting(blob: str) -> str:
    return first_match(SETTING_PATTERNS, blob, "in_vitro")


def classify_delivery(blob: str, kind: str) -> str:
    stated = first_match(DELIVERY_PATTERNS, blob, "unknown")
    if re.search(r"intestinal exposure unknown|oral gut efficacy unresolved|delivery (?:evidence )?"
                 r"(?:is )?(?:none|unmeasured|untested|not measured)", blob.lower()):
        stated = "unknown"
    return DEL.delivery_for_kind(kind, stated=None if stated == "unknown" else stated)


def classify_direction(blob: str) -> str:
    return first_match(DIRECTION_PATTERNS, blob, "unclear")


def direction_near(name: str, blob: str, *, window: int = 90) -> str | None:
    """Which way this one organism went, read from the words nearest its name.

    A cascade row describes several organisms at once - "R. gnavus relative
    enrichment; B. uniformis, B. thetaiotaomicron and A. muciniphila
    decline" - and one direction for the whole row is wrong for at least one
    of them. So each name gets its own window, and within that window the
    direction word *closest to the name* wins rather than whichever pattern
    happens to be checked first: "enrichment" sits beside R. gnavus and
    "decline" sits beside the three organisms it is actually about.
    """
    genus, _, epithet = name.partition(" ")
    patterns = [re.escape(name)]
    if epithet and genus:
        patterns.append(rf"\b{re.escape(genus[0])}\.\s*{re.escape(epithet)}")
    best: tuple[int, str] | None = None
    for pattern in patterns:
        for m in re.finditer(pattern, blob, flags=re.IGNORECASE):
            left, right = max(0, m.start() - window), m.end() + window
            for direction, expr in DIRECTION_PATTERNS:
                for hit in re.finditer(expr, blob[left:right], flags=re.IGNORECASE):
                    at = left + hit.start()
                    distance = 0 if m.start() <= at <= m.end() else min(
                        abs(at - m.end()), abs(m.start() - at))
                    if best is None or distance < best[0]:
                        best = (distance, direction)
    return best[1] if best else None


def direction_of(cells: dict[str, str], whole: str) -> str:
    """Which way the material moved the target, read from the nearest cell.

    Reading the whole row at once inverted findings: a boundary column that
    says "does not establish eradication" sits beside a result that says the
    growth increased, and whichever pattern happened to match first won. So
    the cell that states the result is asked first, then the target, and the
    rest of the row only when neither says anything. A row that still says
    nothing stays `unclear`, and an unclear row is never offered as a way to
    fix anything.
    """
    for cell in (cells.get("effect", ""), cells.get("target", ""), cells.get("material", "")):
        if not cell:
            continue
        found = classify_direction(cell)
        if found != "unclear":
            return found
    return classify_direction(whole)


# --------------------------------------------------------------------------- #
# targets
# --------------------------------------------------------------------------- #

_BINOMIAL = re.compile(r"\b([A-Z][a-z]{2,})\s+([a-z][a-z-]{2,})\b")
_ABBREV = re.compile(r"\b([A-Z])\.\s*([a-z][a-z-]{2,})\b")

#: Abbreviated species names are everywhere in this literature, and expanding
#: them by guessing the genus from the initial alone invents organisms:
#: "C. difficile" became Candida difficile, Collinsella difficile and
#: Citrobacter difficile at once. The epithet is what identifies the species,
#: so an abbreviation is expanded only when its epithet is known and its
#: initial agrees, or when the row itself names the genus in full somewhere.
EPITHET_GENUS: dict[str, str] = {
    "difficile": "Clostridioides", "coli": "Escherichia", "gnavus": "Ruminococcus",
    "pneumoniae": "Klebsiella", "fragilis": "Bacteroides", "muciniphila": "Akkermansia",
    "prausnitzii": "Faecalibacterium", "aureus": "Staphylococcus", "pylori": "Helicobacter",
    "albicans": "Candida", "boulardii": "Saccharomyces", "cerevisiae": "Saccharomyces",
    "nucleatum": "Fusobacterium", "faecalis": "Enterococcus", "faecium": "Enterococcus",
    "lenta": "Eggerthella", "vulgatus": "Bacteroides", "thetaiotaomicron": "Bacteroides",
    "uniformis": "Bacteroides", "distasonis": "Parabacteroides", "aeruginosa": "Pseudomonas",
    "subtilis": "Bacillus", "butyricum": "Clostridium", "perfringens": "Clostridium",
    "longum": "Bifidobacterium", "breve": "Bifidobacterium", "bifidum": "Bifidobacterium",
    "adolescentis": "Bifidobacterium", "infantis": "Bifidobacterium", "animalis": "Bifidobacterium",
    "lactis": "Bifidobacterium", "reuteri": "Limosilactobacillus", "rhamnosus": "Lacticaseibacillus",
    "plantarum": "Lactiplantibacillus", "acidophilus": "Lactobacillus", "casei": "Lacticaseibacillus",
    "hadrus": "Anaerostipes", "rectale": "Agathobacter", "eligens": "Lachnospira",
    "bromii": "Ruminococcus", "wadsworthia": "Bilophila", "smithii": "Methanobrevibacter",
    "scindens": "Clostridium", "innocuum": "Clostridium", "hiranonis": "Peptacetobacter",
    "copri": "Prevotella", "typhimurium": "Salmonella", "enterica": "Salmonella",
    "jejuni": "Campylobacter", "monocytogenes": "Listeria", "mirabilis": "Proteus",
    "freundii": "Citrobacter", "cloacae": "Enterobacter", "aerofaciens": "Collinsella",
    "parvula": "Veillonella", "intestinalis": "Bacteroides", "hallii": "Anaerobutyricum",
}

#: Genera the atlas abbreviates, by initial, used only to check that an
#: expansion is consistent with the letter the source actually wrote.
INITIALS: dict[str, tuple[str, ...]] = {
    "A": ("Akkermansia", "Anaerostipes", "Alistipes", "Aspergillus", "Anaerotignum", "Agathobacter",
          "Anaerobutyricum"),
    "B": ("Bacteroides", "Bifidobacterium", "Blautia", "Bacillus", "Bilophila"),
    "C": ("Clostridium", "Clostridioides", "Candida", "Collinsella", "Citrobacter", "Campylobacter"),
    "D": ("Desulfovibrio", "Dorea"),
    "E": ("Escherichia", "Enterococcus", "Eubacterium", "Enterobacter", "Eggerthella"),
    "F": ("Faecalibacterium", "Fusobacterium", "Flavonifractor"),
    "H": ("Helicobacter",),
    "K": ("Klebsiella",),
    "L": ("Lactobacillus", "Lacticaseibacillus", "Limosilactobacillus", "Lactiplantibacillus",
          "Lachnospira", "Listeria"),
    "M": ("Mediterraneibacter", "Methanobrevibacter"),
    "P": ("Prevotella", "Parabacteroides", "Phocaeicola", "Pseudomonas", "Proteus", "Parasutterella",
          "Peptacetobacter"),
    "R": ("Ruminococcus", "Roseburia"),
    "S": ("Streptococcus", "Staphylococcus", "Salmonella", "Saccharomyces", "Shigella"),
    "V": ("Veillonella", "Vibrio"),
}

#: Words that look like a binomial but are not organisms.
_NOT_A_TAXON = {
    "in vitro", "ex vivo", "human milk", "primary source", "randomized trial", "high fat",
    "short chain", "bile acid", "fatty acid", "gut microbiome", "small intestine", "large intestine",
    "clinical trial", "case report", "whole grain", "resistant starch", "essential oil",
}

#: Report panels, by the words a source would use for them.
#:
#: Derived from the panel definitions themselves rather than written out by
#: hand. A hand-kept list covered 29 of the 46 panels, so seventeen readings
#: - agmatine, B6, glutathione, spermidine, thiamine and the rest - could
#: never be matched to a single piece of evidence no matter what the research
#: said. Each panel contributes its own name, its metabolite, its pathway and
#: the gene symbols of its targets, which is what a paper actually calls these
#: things.
PANELS_DIR = REPO / "panels"

#: Words too generic to match on: they would attach half the atlas to a panel.
_PANEL_STOPWORDS = {
    "formation", "production", "capacity", "synthesis", "degradation", "metabolism",
    "pathway", "routes", "route", "activity", "gut", "microbial", "bacterial", "total",
    "and", "the", "of", "acid", "acids", "vitamin", "amino", "other", "related",
}


def _panel_terms() -> dict[str, str]:
    """A regex per panel, built from how the panel itself is defined."""
    out: dict[str, str] = {}
    for path in sorted(PANELS_DIR.glob("*.yaml")):
        if path.name.startswith("_"):
            continue
        try:
            spec = yaml.safe_load(path.read_text()) or {}
        except (OSError, yaml.YAMLError):
            continue
        panel = path.stem
        words: set[str] = set()
        for value in (spec.get("name"), spec.get("metabolite"), spec.get("pathway")):
            for part in re.split(r"[^A-Za-z0-9'-]+", str(value or "")):
                part = part.strip("'-")
                if len(part) > 3 and part.lower() not in _PANEL_STOPWORDS:
                    words.add(part.lower())
        for target in (spec.get("targets") or ()):
            gene = str((target or {}).get("gene") or "").strip()
            if len(gene) >= 3:
                words.add(gene.lower())
        if panel not in _PANEL_STOPWORDS and len(panel) > 2:
            words.add(panel.lower())
        if words:
            out[panel] = r"\b(?:" + "|".join(re.escape(w) for w in sorted(words)) + ")"
    return out


#: The hand-written patterns stay as overrides where a panel needs a synonym
#: its own definition does not contain - "p-cresol" for pcresol, "TMAO" for
#: cutc - and are merged over the derived ones.
PANEL_OVERRIDES: dict[str, str] = {
    "butyrate": r"\b(?:butyrate|butyric|butyrogen)",
    "propionate": r"\b(?:propionate|propionic)",
    "bsh": r"\b(?:bile salt hydrolase|\bbsh\b|deconjugat)",
    "bai": r"\b(?:secondary bile acid|7-?alpha-?dehydroxyl|\bbai\b|deoxychol|lithochol)",
    "h2s": r"\b(?:hydrogen sulf|sulphide|sulfide|\bh2s\b|sulfate-?reduc|sulphate-?reduc)",
    "methane": r"\b(?:methane|methanogen|archaea)",
    "bglucuronidase": r"\b(?:beta-?glucuronidase|\bgus\b|glucuronide)",
    "cutc": r"\b(?:choline|\bcutc\b|\btmao\b|trimethylamine)",
    "pcresol": r"\b(?:p-?cresol|cresyl|cresol)",
    "indole": r"\b(?:indole)",
    "ipa": r"\b(?:indole-?3-?propion|\bipa\b)",
    "tryptamine": r"\b(?:tryptamine|tryptophan)",
    "histamine": r"\b(?:histamine|histidine decarboxyl)",
    "gaba": r"\b(?:\bgaba\b|gamma-?aminobutyr|glutamate decarboxyl)",
    "dopamine": r"\b(?:dopamine|levodopa|tyrosine decarboxyl)",
    "b12": r"\b(?:\bb12\b|cobalamin|cbi\b)",
    "folate": r"\b(?:folate|folic acid)",
    "riboflavin": r"\b(?:riboflavin|vitamin b2)",
    "biotin": r"\b(?:biotin|vitamin b7)",
    "k2": r"\b(?:vitamin k2|menaquinone)",
    "oxalate": r"\b(?:oxalate|oxalobacter)",
    "urease": r"\b(?:urease|ammonia)",
    "equol": r"\b(?:equol|isoflavone|daidzein|genistein)",
    "urda": r"\b(?:urolithin|ellagitannin|ellagic)",
    "mucinlps": r"\b(?:mucin|mucus layer|lipopolysacchar|\blps\b|lipid a\b)",
    "crc_virulence": r"\b(?:colibactin|\bpks\b|fragilysin|\bbft\b|virulence)",
    "carbohydrates": r"\b(?:carbohydrate-?active|\bcazy|glycoside hydrolase|polysacchar)",
    "fermentation": r"\b(?:fermentation|short-?chain fatty acid|\bscfa\b)",
    "bcaa": r"\b(?:branched-?chain|\bbcaa\b|valine|leucine|isoleucine)",
    "glutathione": r"\b(?:glutathione)",
    "spermidine": r"\b(?:spermidine|polyamine)",
    "putrescine": r"\b(?:putrescine|polyamine)",
    "ethanol": r"\b(?:ethanol|alcohol)",
}

PANEL_WORDS: dict[str, str] = {**_panel_terms(), **PANEL_OVERRIDES}

#: Report-level concerns that are not organisms or panels.
CONCERN_WORDS: dict[str, str] = {
    # Whole-gut subjects a report raises that are not one organism and not
    # one pathway. Seven of these left 120 atlas rows bound to nothing at
    # all - real materials with real papers that no finding could ever reach.
    "biofilm": r"\bbiofilm|\beps\b|extracellular polymeric",
    "barrier": r"\bbarrier|tight junction|permeab|leaky gut|mucosal (?:support|heal)|"
               r"epithelial integrity|zonulin|occludin",
    "mucus": r"\bmucus|mucin|goblet cell",
    "inflammation": r"\binflammat|calprotectin|\btnf\b|il-?6|il-?10|cytokine|crp\b|"
                    r"nf-?kappa|colitis",
    "pathogen_burden": r"\bpathogen|decolon|eradicat|infection|virulence|toxin",
    "diversity": r"\bdiversity|richness|evenness|shannon|dysbio",
    "transit": r"\bconstipat|diarrh|stool (?:form|consistency|frequency)|transit|motility|bloat",
    "antibiotic_resistance": r"\bresistance gene|\bamr\b|resistome|multidrug-?resistant|\besbl\b",
    "immune": r"\bimmune|iga\b|treg|t-?cell|vaccine response|allerg",
    "metabolic": r"\binsulin|glycaem|glycem|glucose tolerance|hba1c|adiposity|body weight|"
                 r"lipid|cholesterol|triglycerid",
    "gut_brain": r"\bcognit|mood|depress|anxiet|fatigue|brain|sleep",
    "colonisation_resistance": r"\bcolonization resistance|colonisation resistance|"
                               r"competitive exclusion|niche exclusion|engraft",
    "bile_acids": r"\bbile acid|bile salt|deconjugat|cholic|chenodeoxychol",
    "short_chain_fatty_acids": r"\bshort-?chain fatty acid|\bscfa\b|butyrate|propionate|acetate",
    "fungal": r"\bcandida|fungal|mycobiome|yeast|aspergillus|saccharomyces",
}


def taxa_in(blob: str) -> list[str]:
    """Every organism a cell names, normalised to `Genus species`.

    An abbreviation expands only to the genus its epithet actually belongs
    to, or to a genus the same text spells out in full. Anything else stays
    unexpanded rather than becoming an organism nobody studied.
    """
    out: list[str] = []
    spelled = {m.group(1) for m in _BINOMIAL.finditer(blob)}
    for m in _BINOMIAL.finditer(blob):
        name = f"{m.group(1)} {m.group(2)}"
        if name.lower() in _NOT_A_TAXON or m.group(1) in {"The", "This", "Not", "No", "In", "Human", "Both"}:
            continue
        out.append(name)
    for m in _ABBREV.finditer(blob):
        initial, epithet = m.group(1), m.group(2)
        known = EPITHET_GENUS.get(epithet)
        if known and known.startswith(initial):
            out.append(f"{known} {epithet}")
            continue
        here = [g for g in INITIALS.get(initial, ()) if g in spelled]
        if len(here) == 1:
            out.append(f"{here[0]} {epithet}")
    seen: list[str] = []
    for n in out:
        if n not in seen:
            seen.append(n)
    return seen


def genera_in(blob: str) -> list[str]:
    """Genus-level mentions: `Bifidobacterium`, `Lactobacillus species`."""
    out = []
    for m in re.finditer(r"\b([A-Z][a-z]{4,})\s*(?:spp?\.|species|genus)\b", blob):
        out.append(m.group(1))
    seen: list[str] = []
    for n in out:
        if n not in seen:
            seen.append(n)
    return seen


def matches_in(words: dict[str, str], blob: str) -> list[str]:
    low = blob.lower()
    return [key for key, pattern in words.items() if re.search(pattern, low)]


# --------------------------------------------------------------------------- #
# the import
# --------------------------------------------------------------------------- #


def rows_of(cur: sqlite3.Cursor, table_id: str) -> list[dict[str, Any]]:
    out = []
    for legacy, fj, line in cur.execute(
        "select legacy_id, fields_json, source_line from records where table_id=? order by row_index",
        (table_id,),
    ):
        out.append({"legacy_id": legacy, "fields": json.loads(fj), "source_line": line})
    return out


def option_from(row: dict[str, Any], table_id: str, roles: dict[str, str]) -> dict[str, Any] | None:
    fields = row["fields"]
    cells = {role: text_of(fields.get(col, "")) for col, role in roles.items() if role != "id"}
    material = clean_material(cells.get("material", ""))
    if not material or len(material) < 3:
        return None
    whole = " ".join(text_of(v) for v in fields.values())
    target_blob = " ".join(x for x in (cells.get("target", ""), cells.get("effect", "")) if x)
    if not target_blob:
        target_blob = whole

    taxa = taxa_in(target_blob)[:8]
    kind = classify_kind(material, whole)
    setting = classify_setting(" ".join(x for x in (cells.get("setting", ""), whole) if x))
    route = classify_delivery(" ".join(x for x in (cells.get("delivery", ""), cells.get("setting", ""), whole) if x), kind)
    reach = DEL.assess(setting, route)
    if not reach.reportable:
        return None

    source_blob = html.unescape(str(fields.get(next((c for c, r in roles.items() if r == "source"), ""), "")))
    urls = urls_of(source_blob) or urls_of(whole)
    effect = cells.get("effect") or cells.get("target") or cells.get("setting") or ""
    option = {
        "option_id": "ATLAS-" + str(row["legacy_id"] or f"{table_id}-{row['source_line']}"),
        "source_table": table_id,
        "material": material[:240],
        "kind": kind,
        "setting": setting,
        "delivery": route,
        "direction": direction_of(cells, whole),
        "effect": (effect or "")[:400],
        "boundary": (cells.get("boundary") or "")[:400],
        "taxa": taxa,
        "taxon_direction": {
            t: d for t in taxa
            if (d := direction_near(t, " ".join(x for x in (cells.get("effect", ""),
                                                            cells.get("target", ""),
                                                            cells.get("material", "")) if x))) is not None
        },
        "genera": genera_in(target_blob)[:6],
        "panels": matches_in(PANEL_WORDS, target_blob),
        "concerns": matches_in(CONCERN_WORDS, target_blob),
    }
    if urls:
        option["url"] = urls[0]
    if (pmid := pmid_of(source_blob)):
        option["pmid"] = pmid
    if (doi := doi_of(source_blob)):
        option["doi"] = doi
    return option


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--zip", type=Path, default=DEFAULT_ZIP)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--report", action="store_true", help="print what was imported and stop")
    args = ap.parse_args()

    if not args.zip.is_file():
        sys.exit(f"{args.zip} not found")
    with tempfile.TemporaryDirectory() as tmp:
        with zipfile.ZipFile(args.zip) as zf:
            name = next(n for n in zf.namelist() if n.endswith(".sqlite"))
            db = Path(zf.extract(name, tmp))
        con = sqlite3.connect(f"file:{db.resolve().as_uri()[7:]}?mode=ro", uri=True)
        head, body = con.cursor(), con.cursor()

        options: list[dict[str, Any]] = []
        per_table: dict[str, int] = {}
        for table_id in (*EVIDENCE_TABLES, *CASCADE_TABLES):
            found = head.execute(
                "select headers_json from data_tables where table_id=?", (table_id,)
            ).fetchone()
            if not found:
                continue
            headers = json.loads(found[0])
            roles: dict[str, str] = {}
            for col in headers:
                role = role_of(col)
                if role and role not in set(roles.values()):
                    roles[col] = role
            if "material" not in roles.values():
                # the first non-ID column is the material in every atlas table
                for col in headers:
                    if role_of(col) != "id":
                        roles[col] = "material"
                        break
            kept = 0
            for row in rows_of(body, table_id):
                option = option_from(row, table_id, roles)
                if option is not None:
                    options.append(option)
                    kept += 1
            per_table[table_id] = kept
        con.close()

    settings = collections.Counter(o["setting"] for o in options)
    kinds = collections.Counter(o["kind"] for o in options)
    deliveries = collections.Counter(o["delivery"] for o in options)
    directions = collections.Counter(o["direction"] for o in options)
    bound = sum(1 for o in options if o["taxa"] or o["genera"] or o["panels"] or o["concerns"])

    print(f"imported {len(options)} options from {len(per_table)} tables")
    print(f"  bound to at least one target: {bound} ({100 * bound // max(len(options), 1)}%)")
    print(f"  with a citation: {sum(1 for o in options if o.get('url') or o.get('pmid') or o.get('doi'))}")
    print(f"  settings:   {dict(settings.most_common())}")
    print(f"  kinds:      {dict(kinds.most_common())}")
    print(f"  delivery:   {dict(deliveries.most_common())}")
    print(f"  direction:  {dict(directions.most_common())}")
    if args.report:
        for o in options[:6]:
            print("   ", json.dumps(o)[:220])
        return 0

    import yaml  # noqa: PLC0415

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        "# Generated by scripts/import_atlas_evidence.py from\n"
        "# specs/research/OPENBIOTA_INTERVENTION_DATA.zip. Do not edit by hand.\n"
        "#\n"
        "# Every study design is here - human, animal and laboratory culture alike.\n"
        "# `setting` and `delivery` are what the report uses to decide what has to\n"
        "# be said beside a finding; see openbiota/delivery.py.\n"
        + yaml.safe_dump({"atlas_options": options}, sort_keys=False, allow_unicode=True, width=100)
    )
    print(f"wrote {args.out.relative_to(REPO)} ({args.out.stat().st_size / 1024:.0f} KB)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
