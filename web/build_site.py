#!/usr/bin/env python3
"""Build the website in web/openbiota.com: documentation, icons, sitemap, robots, llms.txt.

    web/build_site.py              everything
    web/build_site.py --no-docs    skip the mkdocs build (icons, sitemap, robots, head tags only)
    web/build_site.py --check      report what is missing or stale, change nothing

Everything the site needs lives under web/. The one outside input is the
documentation source (docs/*.md + mkdocs.yml), found through DOCS_SOURCE
(default: the parent of web/, i.e. this repository; set it to another
checkout when the site moves to its own repository).

The landing pages (`index.html`, `research.html`) are hand-written; the
documentation under `docs/` is generated from `docs/*.md` by mkdocs. This
script ties them into one publishable site:

* one `sitemap.xml` at the root that lists every page - the landing pages
  and every documentation page - with the date each source last changed;
* `robots.txt` pointing at that sitemap;
* icons in every size a browser or phone asks for, rasterised from the
  vector mark, plus `site.webmanifest`;
* the `<head>` of each landing page carrying a canonical URL, the icon
  links and Open Graph / Twitter cards.

mkdocs writes its own sitemap inside `docs/`; it is removed after the build
so the root sitemap is the only one.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

WEB = Path(__file__).resolve().parent
SITE = WEB / "openbiota.com"
#: The one input that lives outside web/: the documentation sources and the
#: mkdocs config that renders them into openbiota.com/docs/. Point
#: DOCS_SOURCE at another checkout when the site has its own repository.
DOCS_SOURCE = Path(os.environ.get("DOCS_SOURCE", WEB.parent)).resolve()
BASE_URL = "https://openbiota.com"
MARK = SITE / "img" / "openbiota-mark.svg"
OG_IMAGE = "img/openbiota-og.jpg"
ICON_SIZES = {"favicon-16x16.png": 16, "favicon-32x32.png": 32, "apple-touch-icon.png": 180,
              "android-chrome-192x192.png": 192, "android-chrome-512x512.png": 512}
THEME_COLOUR = "#092b22"


def _git_date(path: Path) -> str:
    """Date of the last commit touching `path`, else the file's mtime."""
    try:
        out = subprocess.run(["git", "log", "-1", "--format=%cs", "--", str(path)], cwd=path.parent,
                             capture_output=True, text=True, check=False).stdout.strip()
        if out:
            return out
    except OSError:
        pass
    return dt.date.fromtimestamp(path.stat().st_mtime).isoformat() if path.exists() else dt.date.today().isoformat()


def _mkdocs() -> list[str]:
    """The mkdocs to run: the docs checkout's venv if it has one, else whatever is on PATH."""
    venv_mkdocs = DOCS_SOURCE / ".venv" / "bin" / "mkdocs"
    if venv_mkdocs.is_file():
        return [str(venv_mkdocs)]
    if shutil.which("mkdocs"):
        return ["mkdocs"]
    return [sys.executable, "-m", "mkdocs"]


def build_docs() -> None:
    """Render DOCS_SOURCE/docs/*.md into the site's docs/ with the checkout's mkdocs.yml (strict)."""
    config = DOCS_SOURCE / "mkdocs.yml"
    if not config.is_file():
        sys.exit(f"no mkdocs.yml at {DOCS_SOURCE}; set DOCS_SOURCE to the checkout that holds docs/ and mkdocs.yml")
    subprocess.run([*_mkdocs(), "build", "--strict", "--config-file", str(config), "--site-dir", str(SITE / "docs")],
                   cwd=DOCS_SOURCE, check=True)
    for stray in ("sitemap.xml", "sitemap.xml.gz"):
        (SITE / "docs" / stray).unlink(missing_ok=True)
    # a markdown twin beside every page, for agents (llms.txt points at these)
    for md in sorted((DOCS_SOURCE / "docs").glob("*.md")):
        html = SITE / "docs" / (md.stem.lower() + ".html")
        if html.is_file():
            shutil.copyfile(md, html.with_suffix(".html.md"))


def build_icons(*, check: bool = False) -> list[str]:
    missing = [name for name in (*ICON_SIZES, "favicon.ico", "site.webmanifest") if not (SITE / name).is_file()]
    if check or not missing:
        return missing
    rsvg = shutil.which("rsvg-convert")
    magick = shutil.which("magick") or shutil.which("convert")
    if rsvg is None:
        raise SystemExit("rsvg-convert is needed to rasterise the icons (brew install librsvg)")
    for name, px in ICON_SIZES.items():
        subprocess.run([rsvg, "-w", str(px), "-h", str(px), "-o", str(SITE / name), str(MARK)], check=True)
    if magick:
        subprocess.run([magick, str(SITE / "favicon-16x16.png"), str(SITE / "favicon-32x32.png"),
                        str(SITE / "favicon.ico")], check=True)
    manifest = {
        "name": "OpenBiota", "short_name": "OpenBiota",
        "description": "Free, open gut microbiome reports from shotgun metagenomic reads.",
        "icons": [{"src": "/android-chrome-192x192.png", "sizes": "192x192", "type": "image/png"},
                  {"src": "/android-chrome-512x512.png", "sizes": "512x512", "type": "image/png"}],
        "theme_color": THEME_COLOUR, "background_color": "#ffffff", "display": "standalone", "start_url": "/",
    }
    (SITE / "site.webmanifest").write_text(json.dumps(manifest, indent=2) + "\n")
    return []


HEAD_TAGS = """<link rel="canonical" href="{url}">
<link rel="icon" href="/favicon.ico" sizes="32x32">
<link rel="icon" type="image/png" sizes="32x32" href="/favicon-32x32.png">
<link rel="icon" type="image/png" sizes="16x16" href="/favicon-16x16.png">
<link rel="apple-touch-icon" sizes="180x180" href="/apple-touch-icon.png">
<link rel="manifest" href="/site.webmanifest">
<meta property="og:type" content="website">
<meta property="og:site_name" content="OpenBiota">
<meta property="og:url" content="{url}">
<meta property="og:title" content="{title}">
<meta property="og:description" content="{description}">
<meta property="og:image" content="{base}/{og_image}">
<meta property="og:image:type" content="image/jpeg">
<meta property="og:image:width" content="1200">
<meta property="og:image:height" content="630">
<meta name="twitter:card" content="summary_large_image">
<meta name="twitter:title" content="{title}">
<meta name="twitter:description" content="{description}">
<meta name="twitter:image" content="{base}/{og_image}">
<meta name="robots" content="index,follow">"""


def _attr(html: str, pattern: str) -> str:
    m = re.search(pattern, html, re.S)
    return (m.group(1) if m else "").strip()


def add_head_tags(page: Path, *, check: bool = False) -> bool:
    """Canonical, icons and social cards in one landing page's head; idempotent."""
    html = page.read_text(encoding="utf-8")
    if 'rel="canonical"' in html:
        return False
    if check:
        return True
    url = f"{BASE_URL}/" if page.name == "index.html" else f"{BASE_URL}/{page.name}"
    title = _attr(html, r"<title>(.*?)</title>")
    description = _attr(html, r'<meta name="description" content="(.*?)"')
    tags = HEAD_TAGS.format(url=url, title=title.replace('"', "&quot;"), description=description,
                            base=BASE_URL, og_image=OG_IMAGE)
    # the vector icon stays; the rasters are added beside it
    anchor = '<link rel="icon" type="image/svg+xml" href="img/openbiota-mark.svg">'
    if anchor in html:
        html = html.replace(anchor, anchor + tags, 1)
    else:
        html = html.replace("</head>", tags + "</head>", 1)
    page.write_text(html, encoding="utf-8")
    return True


def build_sitemap() -> int:
    entries: list[tuple[str, str]] = []
    for name in ("index.html", "research.html"):
        page = SITE / name
        if page.is_file():
            loc = f"{BASE_URL}/" if name == "index.html" else f"{BASE_URL}/{name}"
            entries.append((loc, _git_date(page)))
    docs_out = SITE / "docs"
    for page in sorted(docs_out.glob("*.html")):
        if page.name == "404.html":
            continue
        source = DOCS_SOURCE / "docs" / (page.stem.upper() + ".md")
        if not source.is_file():
            source = DOCS_SOURCE / "docs" / (page.stem + ".md")
        entries.append((f"{BASE_URL}/docs/{page.name}", _git_date(source if source.is_file() else page)))
    lines = ['<?xml version="1.0" encoding="UTF-8"?>',
             '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">']
    for loc, lastmod in entries:
        priority = "1.0" if loc.endswith("openbiota.com/") else ("0.8" if "/docs/" not in loc else "0.6")
        lines += ["  <url>", f"    <loc>{loc}</loc>", f"    <lastmod>{lastmod}</lastmod>",
                  f"    <priority>{priority}</priority>", "  </url>"]
    lines.append("</urlset>")
    (SITE / "sitemap.xml").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (SITE / "robots.txt").write_text(f"User-agent: *\nAllow: /\n\nSitemap: {BASE_URL}/sitemap.xml\n# Agents: a curated map of this site is at {BASE_URL}/llms.txt\n", encoding="utf-8")
    return len(entries)


#: One line per documentation page for llms.txt, keyed by source stem.
_DOC_NOTES = {
    "index": "what the software does, in one page",
    "QUICKSTART": "from FASTQ files to a report in a few commands",
    "INSTALL": "tools, reference data (sizes, sources), hardware and runtime",
    "USAGE": "every command-line flag of `openbiota run` and the other subcommands",
    "OUTPUT": "the report section by section, and every field of results.json",
    "METHOD": "how each reading is computed",
    "PANELS": "the 25 metabolite-pathway panels and their gene families",
    "PROFILES": "the disease-pattern profiles, their sources and how they are scored",
    "FMT": "donor matching for faecal microbiota transplantation",
    "EXPANDED_DETECTION": "the nine detection lanes, competitive confirmation, one share per organism, levels against the reference",
    "MEASURED_CAPABILITY": "recall, precision and false-positive rates measured on simulated communities with known truth",
    "VALIDATION": "cross-platform concordance and the benchmark design",
    "TESTING": "the test suite and what it guarantees",
    "ATTRIBUTION": "every reference database and tool, with licences",
}


def build_llms_txt() -> None:
    """/llms.txt and /docs/llms.txt per https://llmstxt.org/ - a short, curated
    map of the site for agents, pointing at the markdown twins of the docs."""
    docs_dir = SITE / "docs"
    pages = []
    for stem, note in _DOC_NOTES.items():
        md = docs_dir / (stem.lower() + ".html.md")
        if md.is_file():
            title = next((ln.lstrip("# ").strip() for ln in md.read_text(encoding="utf-8").splitlines()
                          if ln.startswith("# ")), stem.replace("_", " ").title())
            pages.append(f"- [{title}]({BASE_URL}/docs/{md.name}): {note}")
    head = [
        "# OpenBiota Gut Health Test",
        "",
        "> Open-source software that turns the paired-end FASTQ files of a shotgun stool metagenome into a 300+ page "
        "gut-health report: nine detection methods over the current reference catalogues, one reconciled share per "
        "organism, levels against the reference adults who carry each organism, a 485-target pathogen screen, 25 "
        "metabolite-pathway capacities, disease-pattern resemblance, and graded intervention evidence. Every number "
        "in the PDF is also in results.json. Research use only; not a diagnostic test.",
        "",
        "The software is source-available under the PolyForm Noncommercial licence at "
        "https://github.com/openbiota-org/openbiota-gut-health-test. The documentation below is the authoritative "
        "description of what the report contains and how each reading is computed; the markdown files are the same "
        "text as the HTML pages.",
        "",
        "## Documentation",
        "",
        *pages,
        "",
        "## Optional",
        "",
        f"- [Landing page]({BASE_URL}/): what the project is, for a general reader",
        f"- [Research]({BASE_URL}/research.html): the published studies behind each reading",
        "- [Source repository](https://github.com/openbiota-org/openbiota-gut-health-test): code, tests, issues",
        "",
    ]
    text = "\n".join(head)
    (SITE / "llms.txt").write_text(text, encoding="utf-8")
    (docs_dir / "llms.txt").write_text(text, encoding="utf-8")


def add_markdown_alternates() -> int:
    """<link rel="alternate" type="text/markdown"> and rel="describedby" on every docs page."""
    n = 0
    for page in sorted((SITE / "docs").glob("*.html")):
        twin = page.with_suffix(".html.md")
        if not twin.is_file():
            continue
        html = page.read_text(encoding="utf-8")
        if 'type="text/markdown"' in html:
            continue
        tags = (f'<link rel="alternate" type="text/markdown" href="{BASE_URL}/docs/{twin.name}">'
                f'<link rel="describedby" href="{BASE_URL}/docs/llms.txt">')
        html = html.replace("</head>", tags + "</head>", 1)
        page.write_text(html, encoding="utf-8")
        n += 1
    return n


def build_root_404() -> None:
    """CloudFront serves /404.html for a missing object; the docs theme already renders one."""
    src = SITE / "docs" / "404.html"
    if src.is_file():
        html = src.read_text(encoding="utf-8")
        # the docs 404 is built with paths relative to /docs/; make them absolute
        html = html.replace('href="assets/', 'href="/docs/assets/').replace('src="assets/', 'src="/docs/assets/')
        html = html.replace('href="stylesheets/', 'href="/docs/stylesheets/').replace('src="javascripts/', 'src="/docs/javascripts/')
        (SITE / "404.html").write_text(html, encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-docs", action="store_true", help="skip the mkdocs build")
    ap.add_argument("--check", action="store_true", help="report only")
    a = ap.parse_args()
    if a.check:
        missing = build_icons(check=True)
        stale_heads = [p.name for p in (SITE / "index.html", SITE / "research.html") if add_head_tags(p, check=True)]
        for name in ("sitemap.xml", "robots.txt"):
            if not (SITE / name).is_file():
                missing.append(name)
        print("missing:", missing or "nothing", "| heads without canonical/icons:", stale_heads or "none")
        return 1 if (missing or stale_heads) else 0
    if not a.no_docs:
        build_docs()
    build_icons()
    for page in (SITE / "index.html", SITE / "research.html"):
        if page.is_file() and add_head_tags(page):
            print(f"head tags added to {page.name}")
    n = build_sitemap()
    build_llms_txt()
    alternates = add_markdown_alternates()
    build_root_404()
    print(f"sitemap: {n} pages; robots.txt, llms.txt, 404.html written; {alternates} docs pages link their markdown twin; "
          f"icons: {', '.join(ICON_SIZES)} + favicon.ico + site.webmanifest")
    return 0


if __name__ == "__main__":
    sys.exit(main())
