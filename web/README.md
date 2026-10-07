# openbiota.com — the website

Everything that makes and publishes the site lives in this folder, so it can
move to its own repository without touching the rest of the project.

```
web/
  openbiota.com/        the deployed site, exactly as served: index.html, research.html,
                        style.css, script.js, fonts/, img/, docs/ (generated), llms.txt, …
  build_site.py         renders the documentation into openbiota.com/docs/, writes the
                        sitemap, robots.txt, llms.txt, icons, head tags and the 404 page
  deploy_site.sh        builds, checks, uploads to S3 and invalidates CloudFront
  deploy.env            where to deploy (bucket, distribution, zone); not secret
  infra/                one-time AWS setup: Route 53 zone copy, certificate, bucket,
                        CloudFront distribution, DNS aliases
  tools/                img-to-webp.py (image renditions), build_research_data.py
                        (research.html's catalogue from specs/research/)
  img-src/              original full-size images; the site serves the WebP/JPEG renditions
  .venv/                tooling environment (boto3, Pillow); created by `make web-setup`
```

## Editing the landing pages

`index.html`, `research.html`, `style.css`, `script.js` are hand-written and
served as they are — edit them and deploy. No build step applies to them except
the deploy script's cache-busting stamp on `?v=` asset URLs.

Images: put the original in `img-src/`, make the rendition with
`tools/img-to-webp.py` (WebP at display width, quality ~80), reference the
`.webp` from the page with explicit `width`/`height`. The Open Graph image is
`img/openbiota-og.jpg`, 1200×630 JPEG — scrapers (Facebook, LinkedIn, iMessage,
Slack) do not reliably render WebP and reject large files.

Fonts are self-hosted in `fonts/` (DM Sans and Manrope, variable WOFF2, Latin
and Latin-Extended subsets, SIL OFL). The page makes no third-party requests.

## The one input from outside `web/`

The documentation under `/docs/` is rendered from `docs/*.md` with `mkdocs.yml`,
and `research-data.js` from `specs/research/*.md`. Both are found through
`DOCS_SOURCE` (default: the parent of `web/`, i.e. this repository). When the
site has its own repository, set `DOCS_SOURCE` to a checkout of the software
repository:

```
DOCS_SOURCE=~/src/openbiota-gut-health-test web/deploy_site.sh
```

## Deploying

```
web/deploy_site.sh              build, confirm, upload, invalidate
web/deploy_site.sh --dry-run    build and list what would change
web/deploy_site.sh --no-docs    landing-page edits only
```

Needs the AWS CLI and the profile named in `deploy.env`. The script refuses to
deploy if the built site contains a credential-like string or a local sample
identifier. HTML, XML, TXT, JSON and the manifest are uploaded `no-cache`;
CSS, JS, fonts and images for a year (their URLs change when they change).
CloudFront compresses text at the edge, so nothing is pre-compressed.

## Hosting

Route 53 (DNS) → CloudFront (HTTPS, HTTP/2+3, Brotli/gzip, security headers,
all edge locations) → private S3 bucket `openbiota.com` in us-west-1, read
through an Origin Access Control. The certificate covers `openbiota.com` and
`*.openbiota.com`, is DNS-validated and renews itself. `infra/site_setup.py`
created all of it and is idempotent; `infra/route53_zone.py` copied the zone
from the previous DNS host and verifies every record.

Costs: Route 53 $0.50/month; S3 cents; CloudFront within its always-free tier
(1 TB/month). No account other than AWS is involved.
