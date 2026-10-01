# OpenBiota

The complete landing page as plain HTML, CSS, JavaScript, and images.

## Open or deploy

Open `index.html` in your browser, or upload the contents of this folder to
your own repository and any static web host. No installation or build step.
Keep `index.html`, `research.html`, `style.css`, `script.js`, `research.js`,
`research-data.js`, and `img/` together. Relative paths also work when the page
is served from a subfolder.

The design, navigation, expandable sections, sticky buttons, and email modal
work immediately. The original DM Sans and Manrope fonts load from Google
Fonts; offline, the page uses its system-font fallbacks. All images are included.

## The research catalog page

`research.html` lists every published study, dataset, database, and reference
resource used to build the first open-source release, with search, filters, and
per-record provenance notes. The rows come from `research-data.js`, generated
from the research catalogs in `specs/research/` of the project repository —
`OPENBIOTA_PROJECT_RESEARCH_v0.1.0-v0.8.2.md` for the build-specification
history and `OPENBIOTA_PROJECT_RESEARCH_v0.8.3.md` for the current one, merged
so that a source cited by both appears once. `research.js` renders and filters
them in the browser. When a catalog changes, regenerate the data file:

```
python3 tools/build_research_data.py
```

The page itself needs no build step and works from disk, because the data is a
plain script rather than something fetched at runtime.

The home page links to `research.html` from the open science section only —
the “1,300+ published studies” line and “Browse every source” — and those
links open in a new tab so the reader never loses their place. The
institution logos appear twice: as a “Built on published research from” band
closing the report section, which links down to the open science section, and
again, logos only, opening that section. They live in `img/logos/`
(forest/sage, for the light background) and `img/logos/light/` (light green,
for the dark section). They identify research sources only and do not imply
any endorsement or affiliation.

## Section URLs

Each major home-page section has a stable id — `#report`, `#gut-health`,
`#technology`, `#research`, `#story`, `#get-started` — and `script.js` keeps
the address bar on the id of the section under the middle of the viewport as
the reader scrolls (the bare page URL above the first section). The header
links in `index.html`, `research.html`, and the docs template
`docs/overrides/partials/header.html` all point at these ids; change one and
change them all, then rebuild the docs with `make docs-build`.

## Connect the email form

Edit the `OPENBIOTA` settings at the top of `script.js`:

- `repositoryUrl`: the GitHub repository URL (`https://github.com/openbiota-org/openbiota-gut-health-test`).
- `emailEndpoint`: your own form service or server endpoint. Leave blank while
  previewing; no email address is submitted or claimed to be saved.
- `eventsEndpoint`: optional analytics endpoint. Leave blank to disable network
  analytics calls.

The form POSTs JSON with `email`, `interest_kind` (`sample` or `own_report`),
`placement`, `source`, `campaign`, and an empty `website` honeypot field.

The endpoint should validate the email, reject a filled honeypot, save the
request, and return a successful HTTP response with JSON, for example:

```json
{"ok": true}
```

Configure your service to email the sample report and instructions. Only
return `{"ok": true, "delivery": "sent"}` after your service has actually sent
that email; the page then displays “Check your inbox.” Errors should return
an appropriate HTTP error and `{"ok": false, "error": "Please try again."}`.
A third-party endpoint must allow requests from your website through CORS.
Put API keys in your server/service settings, never in these public files.

The sample medical PDF is not included in this website export. Supply the
sample you intend to share through your own delivery service.

## Optional demand measurement

Requests preserve whether a visitor asked for a sample or their own report,
where they clicked, and UTM source/campaign tags. Set `eventsEndpoint` to
receive `landing_view`, `sample_open`, `own_report_open`, and
`report_request_submitted` events as JSON POSTs. Event payloads contain no
email addresses. You can also listen for the browser event
`openbiota:interest` from your own analytics code.

No server, database, dashboard, credentials, or OpenAI hosting configuration
is included or required to display this page. Email capture and delivery are
the only parts that need an external service of your choice.
