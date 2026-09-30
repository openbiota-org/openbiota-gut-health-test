#!/usr/bin/env python3
"""Build the website in web/openbiota.com: documentation, icons, sitemap, robots.

    scripts/build_site.py              everything
    scripts/build_site.py --no-docs    skip the mkdocs build (icons, sitemap, robots, head tags only)
    scripts/build_site.py --check      report what is missing or stale, change nothing

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
import re
import shutil
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SITE = REPO / "web" / "openbiota.com"
BASE_URL = "https://openbiota.com"
MARK = SITE / "img" / "openbiota-mark.svg"
OG_IMAGE = "img/openbiota-ogimg.png"
ICON_SIZES = {"favicon-16x16.png": 16, "favicon-32x32.png": 32, "apple-touch-icon.png": 180,
              "android-chrome-192x192.png": 192, "android-chrome-512x512.png": 512}
THEME_COLOUR = "#092b22"


def _git_date(path: Path) -> str:
    """Date of the last commit touching `path`, else the file's mtime."""
    try:
        out = subprocess.run(["git", "log", "-1", "--format=%cs", "--", str(path)], cwd=REPO,
                             capture_output=True, text=True, check=False).stdout.strip()
        if out:
            return out
    except OSError:
        pass
    return dt.date.fromtimestamp(path.stat().st_mtime).isoformat() if path.exists() else dt.date.today().isoformat()


def build_docs() -> None:
    mkdocs = REPO / ".venv" / "bin" / "mkdocs"
    if not mkdocs.is_file():
        subprocess.run([str(REPO / ".venv" / "bin" / "python"), "-m", "pip", "install", "--quiet", "-e", ".[docs]"],
                       cwd=REPO, check=True)
    subprocess.run([str(mkdocs), "build", "--strict"], cwd=REPO, check=True)
    for stray in ("sitemap.xml", "sitemap.xml.gz"):
        (SITE / "docs" / stray).unlink(missing_ok=True)


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
        source = REPO / "docs" / (page.stem.upper() + ".md")
        if not source.is_file():
            source = REPO / "docs" / (page.stem + ".md")
        entries.append((f"{BASE_URL}/docs/{page.name}", _git_date(source if source.is_file() else page)))
    lines = ['<?xml version="1.0" encoding="UTF-8"?>',
             '<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">']
    for loc, lastmod in entries:
        priority = "1.0" if loc.endswith("openbiota.com/") else ("0.8" if "/docs/" not in loc else "0.6")
        lines += ["  <url>", f"    <loc>{loc}</loc>", f"    <lastmod>{lastmod}</lastmod>",
                  f"    <priority>{priority}</priority>", "  </url>"]
    lines.append("</urlset>")
    (SITE / "sitemap.xml").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (SITE / "robots.txt").write_text(f"User-agent: *\nAllow: /\n\nSitemap: {BASE_URL}/sitemap.xml\n", encoding="utf-8")
    return len(entries)


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
    print(f"sitemap: {n} pages; robots.txt written; icons: {', '.join(ICON_SIZES)} + favicon.ico + site.webmanifest")
    return 0


if __name__ == "__main__":
    sys.exit(main())
