# openbiota.com — the website

Everything that makes and publishes the site is in this folder: the pages, the
documentation theme and config, the build, the deploy, the AWS setup. It is
self-contained, with its own environment, so it stands as a repository of its
own. The one thing it reads from elsewhere is the documentation text, which
lives with the software it documents (see "The one input from outside", below).

```
openbiota.com/        the deployed site, exactly as served: index.html, research.html,
                      style.css, script.js, fonts/, img/, docs/ (generated), llms.txt, …
mkdocs.yml            how the documentation pages are rendered (theme, hooks, nav)
docs-theme/           the docs' look: Material overrides (header, footer, phone bar),
                      the two hooks, stylesheets/openbiota.css
build_site.py         renders the documentation into openbiota.com/docs/, writes the
                      sitemap, robots.txt, llms.txt, icons, head tags and the 404 page;
                      --serve previews the documentation
deploy_site.sh        builds, checks, publishes the edge function, uploads to S3 and
                      invalidates CloudFront
deploy.env            where to deploy (bucket, distribution, zone) and, once the site is
                      its own repository, DOCS_SOURCE; not secret
infra/                one-time AWS setup: Route 53 zone copy, certificate, bucket,
                      CloudFront distribution, the clean-URL function, DNS aliases
tools/                img-to-webp.py (image renditions), build_research_data.py
                      (research.html's catalogue from the software's specs/research/)
img-src/              original full-size images; the site serves the WebP/JPEG renditions
requirements.txt      the tooling: mkdocs, Material, Pillow, boto3
.venv/                tooling environment, ignored:
                      python3 -m venv .venv && .venv/bin/pip install -r requirements.txt
```

## Editing the landing pages

`index.html`, `research.html`, `style.css`, `script.js` are hand-written and
served as they are — edit them and deploy. No build step applies to them except
the deploy script's cache-busting stamp on `?v=` asset URLs.

Images: put the original in `img-src/`, make the renditions with
`tools/img-to-webp.py` (WebP, quality ~80, one file per display width), and
reference them from the page with `srcset`/`sizes` and explicit
`width`/`height` so phones download the small one and nothing shifts. The Open Graph image is
`img/openbiota-og.jpg`, 1200×630 JPEG — scrapers (Facebook, LinkedIn, iMessage,
Slack) do not reliably render WebP and reject large files.

Fonts are self-hosted in `fonts/` (DM Sans and Manrope, variable WOFF2, Latin
and Latin-Extended subsets, SIL OFL). The page makes no third-party requests.

## The one input from outside

The documentation under `/docs/` is rendered from the software repository's
`docs/*.md` (with the screenshots in `docs/images/` and the validation JSON
the pages link to), and `research-data.js` from its `specs/research/*.md`.
Everything about how they are rendered — `mkdocs.yml`, the theme, the hooks,
the stylesheet — is here.

That checkout is `DOCS_SOURCE`. It is read from the environment, else from the
`DOCS_SOURCE` line of `deploy.env`, else it is the folder above this one — the
layout while the site lived inside the software repository. As a repository of
its own, set it once in `deploy.env`; a relative path is relative to this
folder, so with the two repositories side by side:

```
DOCS_SOURCE="../openbiota-gut-health-test"
```

and every script — `build_site.py`, `deploy_site.sh`,
`tools/build_research_data.py` — follows it. `build_site.py --serve` previews
the documentation from the same source at http://127.0.0.1:8000/docs/.

## Deploying

```
./deploy_site.sh              build, confirm, publish the edge function, upload, invalidate
./deploy_site.sh --dry-run    build and list what would change
./deploy_site.sh --no-docs    landing-page edits only
```

This is the only command. It needs the AWS CLI and the profile named in
`deploy.env`, and the tooling environment above (boto3 publishes the edge
function). The script refuses to deploy if the built site contains a
credential-like string or a local sample identifier. HTML, XML, TXT, JSON and the manifest are uploaded `no-cache`;
CSS, JS, fonts and images for a year (their URLs change when they change).
CloudFront compresses text at the edge, so nothing is pre-compressed.

## Hosting

Route 53 (DNS) → CloudFront (HTTPS, HTTP/2+3, Brotli/gzip, security headers,
all edge locations) → private S3 bucket `openbiota.com` in us-west-1, read
through an Origin Access Control. The certificate covers `openbiota.com` and
`*.openbiota.com`, is DNS-validated and renews itself. `infra/site_setup.py`
created all of it and is idempotent; `infra/route53_zone.py` copied the zone
from the previous DNS host and verifies every record.

Clean URLs: the pages link to one another by file name (`index.html`,
`docs/index.html`) so the folder works from disk, but the address bar never
shows `index.html`. A CloudFront Function on viewer requests,
`infra/clean_urls.js`, redirects `/index.html` → `/` and `/docs/index.html` →
`/docs/` (query string kept), adds the slash to `/docs`, and serves `/docs/`
from `docs/index.html` (S3 has no directory index below the root).
`deploy_site.sh` brings it level with the file on every deploy (via
`infra/site_setup.py --clean-urls`: created or updated, run through
CloudFront's test runner, published, attached — a no-op when nothing changed),
so an edit to the function ships like any other change.

Costs: Route 53 $0.50/month; S3 cents; CloudFront within its always-free tier
(1 TB/month). No account other than AWS is involved.
