"""Compile the project research catalogues into the data file behind research.html.

Two catalogues are merged. `OPENBIOTA_PROJECT_RESEARCH_v0.1.0-v0.8.2.md` covers
the build-specification history; `OPENBIOTA_PROJECT_RESEARCH_v0.8.3.md` covers
the v0.8.3 specification and uses its own record IDs, section references and a
much finer type vocabulary. A source cited by both is one record here: the
historical entry keeps its ID, and the newer citation adds its specification
version, section reference, provenance notes, accessions and links. Matching is
by DOI first, then by normalized title.

The verbose per-record type is kept exactly as written, and each record also
gets a coarse family so the page can offer a filter that is short enough to use.

The output is a plain script that assigns one global, not JSON fetched at
runtime, so the page keeps working when `research.html` is opened straight from
disk. Repeated institution, country, type and family names are interned into
lists and referenced by index, and long URL prefixes are shortened;
`research.js` expands both.

Fields the sources leave unestablished are written as empty strings rather than
placeholder words, because the page prints the value as-is.

    python3 tools/build_research_data.py
"""

from __future__ import annotations

import argparse
import collections
import html
import json
import re
import sys
from pathlib import Path

WEB = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(WEB))
from build_site import DOCS_SOURCE  # noqa: E402 - the software checkout, resolved the same way as for the documentation

#: The research catalogues live with the software's specifications.
RESEARCH = DOCS_SOURCE / "specs/research"
HISTORY = RESEARCH / "OPENBIOTA_PROJECT_RESEARCH_v0.1.0-v0.8.2.md"
LATEST = RESEARCH / "OPENBIOTA_PROJECT_RESEARCH_v0.8.3.md"
DESTINATION = WEB / "openbiota.com/research-data.js"

# ---------------------------------------------------------------- catalogues
HISTORY_RECORD_SECTIONS = (
    "Research publications",
    "Datasets, software, databases and supporting resources",
)
HISTORY_NOTES = "Corrections, provenance notes and interpretation limits"
LATEST_RECORD_SECTIONS = (
    "Research publications",
    "Datasets, software, databases and supporting reference resources",
    "Preprints and theses",
    "Product labels and formulation-identity sources",
)
LATEST_NOTES = "Corrections, identity limits and provenance notes"
LATEST_IDENTIFIERS = "Explicit sequence, study and trial identifiers"
LATEST_VERSION = "v0.8.3"

# Source-corpus IDs in the historical catalogue map back to a specification.
SPECIFICATIONS = {
    "C01": "v0.1.0", "C02": "v0.2.0", "C03": "v0.3.0", "C04": "v0.4.0", "C05": "v0.4.1",
    "C06": "v0.4.2", "C07": "v0.5.0", "C08": "v0.6.0", "C09": "v0.7.0", "C10": "v0.8.0",
    "C11": "v0.8.1", "C12": "v0.8.2", "C14": "v0.3.0", "C15": "v0.4.0", "C17": "v0.7.0",
    "C18": "v0.1.0", "C19": "v0.2.0", "C20": "v0.4.0", "C21": "v0.4.0", "C22": "v0.6.0",
    "C23": "v0.7.0", "C24": "v0.7.0",
}
# The remaining source documents are evidence ledgers rather than versions.
DOCUMENTS = {
    "C13": "Disease-evidence ledger",
    "C16": "Donor-scoring specification",
    "C25": "Chronic-urticaria source ledger",
    "C26": "Chronic-urticaria report source",
    "C27": "Claim-source ledger",
    "C28": "Report source ledger",
    "C29": "Long COVID research report",
    "C30": "Disease-profile expansion plan",
    "C31": "Project link collection",
    "C32": "Early citation resolution",
}

# ------------------------------------------------------------------- linking
# Machine-readable endpoints are provenance, not something a reader can use.
API_PATTERNS = (
    "api.crossref.org", "api.datacite.org", "/europepmc/webservices/",
    "/ena/browser/api/", "api.semanticscholar.org", "/api/records",
    "/sviewer/viewer.fcgi", "eutils.ncbi.nlm.nih.gov", "rest.uniprot.org",
)
PRIMARY_LABELS = ("findings", "resource", "publication", "archive", "published version")
URL_PREFIXES = [
    ("1", "https://doi.org/"),
    ("2", "https://europepmc.org/article/MED/"),
    ("3", "https://pubmed.ncbi.nlm.nih.gov/"),
    ("4", "https://pmc.ncbi.nlm.nih.gov/articles/"),
    ("5", "https://www.ncbi.nlm.nih.gov/bioproject/"),
    ("6", "https://github.com/"),
    ("7", "https://www.ebi.ac.uk/ena/browser/view/"),
    ("8", "https://zenodo.org/"),
    ("9", "https://www."),
    ("0", "https://"),
    ("-", "http://"),
]

# --------------------------------------------------------------- type family
# The catalogues describe study designs in their own words — 140+ distinct
# strings. Each record keeps that wording; these families make the filter
# usable. First match wins, so the order is the precedence.
FAMILIES = (
    ("Software, database or method", (
        r"^database", r"^software", r"^method", r"^model data", r"^documentation",
        r"^genome (resource|catalog)", r"^research (database|registry|literature|discovery)",
        r"^reference (database|taxonomy)", r"data infrastructure", r"registry and api")),
    ("Correction or notice", (r"correction", r"retracted", r"publication notice")),
    # \b matters: "evidence synthesis" must not read as a thesis.
    ("Preprint or thesis", (r"preprint", r"\bthesis\b", r"dissertation")),
    ("Product or manufacturer identity", (
        r"label", r"manufacturer", r"developer identity", r"patent")),
    ("Trial registration", (r"trial regist",)),
    ("Guideline or clinical reference", (
        r"guideline", r"consensus", r"regulatory", r"clinical reference",
        r"reporting standard", r"evidence guide", r"evidence summary",
        r"analytical reference report")),
    ("Review or meta-analysis", (r"review", r"meta-analysis", r"evidence synthesis")),
    ("Randomized trial", (r"randomi[sz]ed",)),
    ("Case report", (r"case report", r"case narrative", r"case series")),
    ("Dataset, sequence or code", (
        r"dataset", r"sequence (data|reference)", r"sequencing data", r"supplement",
        r"code", r"research data", r"public sequence")),
    ("Human study", (r"human",)),
    ("Animal study", (r"animal",)),
    ("Laboratory or computational study", (
        r"in vitro", r"ex vivo", r"in vivo", r"biochem", r"computational",
        r"genomic", r"mechanistic", r"mixed")),
    ("Reference or supporting resource", (
        r"supporting resource", r"background article", r"research news")),
)
DEFAULT_FAMILY = "Research publication"

FIELDS = ["id", "name", "year", "description", "contribution", "institutions",
          "countries", "basis", "type", "family", "link", "sources", "notes",
          "specs", "documents", "accession", "isResource"]

BLANK = {"Not established", "Not applicable", "None", "N/A", "", "—", "-"}
LINK_RE = re.compile(r"\[([^\]]+)\]\(<([^>]+)>\)")
YEAR_RE = re.compile(r"\s*\((\d{4})\)\s*$")
ACCESSION_RE = re.compile(r"\s+[—–-]{1,2}\s+((?:[A-Z]{3,5}\d{4,}[,/ ]*)+)$")
TAG_RE = re.compile(r"</?(?:i|b|em|strong|sub|sup)>")
DOI_RE = re.compile(r"doi\.org/(10\.[^\s>)]+)", re.I)


def section(lines: list[str], title: str) -> list[str]:
    """The lines under one `## ` heading."""
    start = next(i + 1 for i, line in enumerate(lines)
                 if line.startswith("## ") and line[3:].strip() == title)
    end = next((j for j in range(start, len(lines)) if lines[j].startswith("## ")), len(lines))
    return lines[start:end]


def rows(lines: list[str], columns: int) -> list[list[str]]:
    """Body rows of a markdown table, with the source's escaped characters decoded."""
    table = []
    for line in lines:
        if not line.startswith("| "):
            continue
        cells = [html.unescape(cell.strip()) for cell in line.strip().strip("|").split("|")]
        if len(cells) != columns or set(cells[0]) <= set("-: "):
            continue
        if cells[0] in ("ID", "Record", "Catalog record", "Accession / identifier"):
            continue
        table.append(cells)
    return table


def shorten(url: str) -> str:
    for code, prefix in URL_PREFIXES:
        if url.startswith(prefix):
            return code + url[len(prefix):]
    return "=" + url


def links(cell: str) -> tuple[str | None, list[str]]:
    """The one link a reader should follow, plus any readable secondary links."""
    primary: str | None = None
    secondary: list[str] = []
    seen: set[str] = set()
    for label, url in LINK_RE.findall(cell):
        url = url.strip()
        if url in seen or any(pattern in url for pattern in API_PATTERNS):
            continue
        seen.add(url)
        if primary is None and label.strip().lower().startswith(PRIMARY_LABELS):
            primary = url
        else:
            secondary.append(url)
    if primary is None and secondary:
        primary = secondary.pop(0)
    return primary, secondary


def family_of(kind: str) -> str:
    lowered = kind.lower()
    for name, patterns in FAMILIES:
        if any(re.search(pattern, lowered) for pattern in patterns):
            return name
    return DEFAULT_FAMILY


def doi_of(cell: str) -> str | None:
    found = DOI_RE.search(cell)
    return found.group(1).lower().rstrip(".").rstrip("/") if found else None


def title_key(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", YEAR_RE.sub("", name).lower())[:90]


def clean(cell: str) -> str:
    return "" if cell in BLANK else TAG_RE.sub("", cell).strip()


def split_values(cell: str) -> list[str]:
    if cell in BLANK:
        return []
    values = []
    for part in (piece.strip() for piece in cell.split(";")):
        if part and part not in BLANK and part not in values:
            values.append(part)
    return values


class Record:
    """One catalogued source, before the fields are interned for the browser."""

    def __init__(self, identifier, name, year, description, contribution,
                 institutions, countries, basis, kind, link, sources, notes,
                 specs, documents, accession, is_resource, doi):
        self.id = identifier
        self.name = name
        self.year = year
        self.description = description
        self.contribution = contribution
        self.institutions = institutions
        self.countries = countries
        self.basis = basis
        self.type = kind
        self.link = link
        self.sources = sources
        self.notes = notes
        self.specs = specs
        self.documents = documents
        self.accession = accession
        self.is_resource = is_resource
        self.doi = doi

    def absorb(self, other: Record) -> None:
        """Fold a later citation of the same source into this record."""
        for note in other.notes:
            if note not in self.notes:
                self.notes.append(note)
        for spec in other.specs:
            if spec not in self.specs:
                self.specs.append(spec)
        for document in other.documents:
            if document not in self.documents:
                self.documents.append(document)
        for url in ([other.link] if other.link else []) + other.sources:
            if url and url != self.link and url not in self.sources:
                self.sources.append(url)
        if other.accession and other.accession not in self.accession:
            self.accession = ", ".join(filter(None, [self.accession, other.accession]))


def version_key(version: str) -> tuple:
    found = re.match(r"v(\d+)\.(\d+)\.(\d+)", version)
    return tuple(int(part) for part in found.groups()) if found else (99, 99, 99)


def read_history(path: Path) -> list[Record]:
    lines = path.read_text(encoding="utf-8").split("\n")
    notes: dict[str, list[str]] = collections.defaultdict(list)
    for record_id, note in rows(section(lines, HISTORY_NOTES), 2):
        note = note.strip()
        if note and note not in notes[record_id]:
            notes[record_id].append(note)

    records = []
    for is_resource, title in enumerate(HISTORY_RECORD_SECTIONS):
        for (record_id, name, description, contribution, institution, country,
             basis, kind, sources, link_cell) in rows(section(lines, title), 10):
            year = ""
            found = YEAR_RE.search(name)
            if found:
                year, name = found.group(1), YEAR_RE.sub("", name)
            accession = ""
            found = ACCESSION_RE.search(name)
            if found:
                accession, name = found.group(1).strip(" ,"), name[: found.start()].strip()
            primary, secondary = links(link_cell)
            source_ids = [source.strip() for source in sources.split(";") if source.strip()]
            records.append(Record(
                record_id, TAG_RE.sub("", name).strip().rstrip(" ."), year,
                clean(description), clean(contribution),
                split_values(institution), split_values(country),
                {"AFFILIATION": "A", "MAINTAINER": "M", "STUDY": "S"}.get(basis, ""),
                kind, primary, secondary, notes.get(record_id, []),
                sorted({SPECIFICATIONS[i] for i in source_ids if i in SPECIFICATIONS},
                       key=version_key),
                sorted({DOCUMENTS[i] for i in source_ids if i in DOCUMENTS}),
                accession, is_resource, doi_of(link_cell)))
    return records


def read_latest(path: Path) -> list[Record]:
    lines = path.read_text(encoding="utf-8").split("\n")
    notes: dict[str, list[str]] = collections.defaultdict(list)
    for record_id, note in rows(section(lines, LATEST_NOTES), 2):
        note = note.strip()
        if note and note not in notes[record_id]:
            notes[record_id].append(note)

    # Accessions are tabled separately and point back at a catalogue record.
    accessions: dict[str, list[str]] = collections.defaultdict(list)
    for (accession, role, data_type, _contribution, record_ref,
         note, _link) in rows(section(lines, LATEST_IDENTIFIERS), 7):
        record_id = record_ref.split()[0].strip(" ,;") if record_ref else ""
        if not re.fullmatch(r"[PRVI]\d{7}", record_id):
            continue
        if accession not in accessions[record_id]:
            accessions[record_id].append(accession)
        detail = " — ".join(filter(None, [f"{accession} ({data_type})" if data_type else accession,
                                          role, clean(note)]))
        if detail and detail not in notes[record_id]:
            notes[record_id].append(detail)

    records = []
    for title in LATEST_RECORD_SECTIONS:
        is_resource = 0 if title == "Research publications" else 1
        for (record_id, name, description, contribution, institution, country,
             basis, kind, sources, link_cell) in rows(section(lines, title), 10):
            year = ""
            found = YEAR_RE.search(name)
            if found:
                year, name = found.group(1), YEAR_RE.sub("", name)
            primary, secondary = links(link_cell)
            # "S083: §15.1 / F30; Explicit identifier" → "v0.8.3 §15.1 / F30".
            documents = []
            for part in sources.split(";"):
                part = part.strip()
                if part.startswith("S083:"):
                    reference = f"{LATEST_VERSION} {part.split(':', 1)[1].strip()}"
                    if reference not in documents:
                        documents.append(reference)
            records.append(Record(
                record_id, TAG_RE.sub("", name).strip().rstrip(" ."), year,
                clean(description), clean(contribution),
                split_values(institution), split_values(country),
                {"AFFILIATION": "A", "MAINTAINER": "M", "STUDY": "S"}.get(basis, ""),
                kind, primary, secondary, notes.get(record_id, []),
                [LATEST_VERSION], documents,
                ", ".join(accessions.get(record_id, [])), is_resource,
                doi_of(link_cell)))
    return records


def merge(history: list[Record], latest: list[Record]) -> tuple[list[Record], int]:
    """Fold a v0.8.3 record into the historical one when both cite the source.

    Only matches against history. Within v0.8.3 a paper, its supplementary data
    and its repository are deliberately separate records that share a DOI, so
    they are left alone. A historical record absorbs at most one repeat; any
    further match stays a record of its own rather than being swallowed.
    """
    by_doi = {record.doi: record for record in reversed(history) if record.doi}
    by_title = {title_key(record.name): record for record in reversed(history)}
    merged = list(history)
    absorbed = 0
    for record in latest:
        existing = by_doi.get(record.doi) if record.doi else None
        if existing is None:
            existing = by_title.get(title_key(record.name))
        if existing is not None:
            existing.absorb(record)
            absorbed += 1
            by_doi.pop(existing.doi, None)
            by_title.pop(title_key(existing.name), None)
            continue
        merged.append(record)
    return merged, absorbed


def build(history_path: Path, latest_path: Path) -> tuple[dict, list[list]]:
    records, absorbed = merge(read_history(history_path), read_latest(latest_path))
    # Published research leads, then datasets; alphabetical inside each group.
    records.sort(key=lambda record: (record.is_resource, record.name.lower()))

    institutions: dict[str, int] = {}
    countries: dict[str, int] = {}
    types: dict[str, int] = {}
    families: dict[str, int] = {}
    documents: dict[str, int] = {}

    def intern(store: dict[str, int], value: str) -> int:
        return store.setdefault(value, len(store))

    payload = []
    for record in records:
        payload.append([
            record.id,
            record.name,
            record.year,
            record.description,
            record.contribution,
            [intern(institutions, name) for name in record.institutions],
            [intern(countries, name) for name in record.countries],
            record.basis,
            intern(types, record.type),
            intern(families, family_of(record.type)),
            shorten(record.link) if record.link else "",
            [shorten(url) for url in record.sources[:4]],
            record.notes,
            sorted(set(record.specs), key=version_key),
            sorted({intern(documents, name) for name in record.documents}),
            record.accession,
            record.is_resource,
        ])

    def names(store: dict[str, int]) -> list[str]:
        return [name for name, _ in sorted(store.items(), key=lambda item: item[1])]

    family_names = names(families)
    tables = {
        "counts": {
            "total": len(payload),
            "publications": sum(1 for row in payload if not row[16]),
            "resources": sum(1 for row in payload if row[16]),
            "institutions": len(institutions),
            # Multinational and International are scopes, not places.
            "countries": len({name for name in countries
                              if name not in ("Multinational", "International")}),
            "linked": sum(1 for row in payload if row[10]),
            "notes": sum(len(row[12]) for row in payload),
            "merged": absorbed,
        },
        "prefixes": dict(URL_PREFIXES),
        "institutions": names(institutions),
        "countries": names(countries),
        "types": names(types),
        "families": family_names,
        "documents": names(documents),
    }
    return tables, payload


def render(tables: dict, records: list[list]) -> str:
    body = ",\n".join(json.dumps(record, ensure_ascii=False, separators=(",", ":"))
                      for record in records)
    return (
        "/* OpenBiota research catalog — generated by tools/build_research_data.py from\n"
        " * specs/research/OPENBIOTA_PROJECT_RESEARCH_v0.1.0-v0.8.2.md and\n"
        " * specs/research/OPENBIOTA_PROJECT_RESEARCH_v0.8.3.md. Do not edit by hand.\n"
        " * Read by research.html. Plain data, no build step or server required.\n"
        " * Record fields: " + ", ".join(FIELDS) + "\n"
        " * institutions/countries/type/family hold indexes into the lists below;\n"
        " * link prefixes are expanded by research.js. */\n"
        "window.OPENBIOTA_RESEARCH = {\n"
        "fields: " + json.dumps(FIELDS) + ",\n"
        "counts: " + json.dumps(tables["counts"], ensure_ascii=False) + ",\n"
        "prefixes: " + json.dumps(tables["prefixes"], ensure_ascii=False) + ",\n"
        "institutions: " + json.dumps(tables["institutions"], ensure_ascii=False) + ",\n"
        "countries: " + json.dumps(tables["countries"], ensure_ascii=False) + ",\n"
        "types: " + json.dumps(tables["types"], ensure_ascii=False) + ",\n"
        "families: " + json.dumps(tables["families"], ensure_ascii=False) + ",\n"
        "documents: " + json.dumps(tables["documents"], ensure_ascii=False) + ",\n"
        "records: [\n" + body + "\n]};\n"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--history", type=Path, default=HISTORY,
                        help="catalogue for v0.1.0–v0.8.2 (default: %(default)s)")
    parser.add_argument("--latest", type=Path, default=LATEST,
                        help="catalogue for the current specification (default: %(default)s)")
    parser.add_argument("--out", type=Path, default=DESTINATION,
                        help="generated data file (default: %(default)s)")
    arguments = parser.parse_args()

    tables, records = build(arguments.history, arguments.latest)
    javascript = render(tables, records)
    arguments.out.write_text(javascript, encoding="utf-8")

    counts = tables["counts"]
    print(f"{arguments.out.relative_to(WEB) if arguments.out.is_relative_to(WEB) else arguments.out}: {counts['total']} records "
          f"({counts['publications']} publications, {counts['resources']} resources), "
          f"{counts['merged']} repeat citations merged, {counts['institutions']} institutions, "
          f"{counts['countries']} countries, {counts['notes']} provenance notes, "
          f"{counts['linked']} linked, {len(javascript.encode()) // 1024} KB")
    missing = [record[0] for record in records if not record[10]]
    if missing:
        print(f"records without a link ({len(missing)}): {', '.join(missing)}")
    print("families: " + ", ".join(
        f"{name} {sum(1 for record in records if record[9] == index)}"
        for index, name in enumerate(tables["families"])))


if __name__ == "__main__":
    main()
