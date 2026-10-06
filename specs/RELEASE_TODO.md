# Release TODO

Everything between here and a public repository, a live site and a release
anyone can install. **[YOU]** needs an account, a decision or a secret;
**[ME]** is code, config or content. Ordered so that each section can start
once the one above it is done; items marked *done* stay for the record.

The operating principle: **the fewest accounts, no recurring cost until there
is real interest, and nothing to run.** Four accounts in total — GitHub (have),
Cloudflare, Brevo, Zenodo — and the only bill is the domain.

---

## 0. Decisions that gate everything else

- **[YOU]** Hosting and DNS: `openbiota.com` on **Cloudflare** (free plan). DNS, static hosting (Pages), inbound mail forwarding (Email Routing), cookieless analytics (Web Analytics), and a backend path later (Workers / Pages Functions on `api.openbiota.com`) in one account. Transfer the registrar to Cloudflare at cost (~$10/yr) or just point the nameservers. *Decided against GitHub Pages: static-only, no headers/redirects/functions, and it was a placeholder I chose for convenience, not a product decision.*
- **[YOU]** Email service for the interest list: **Brevo** — the one service priced by emails *sent*, not contacts *held*: unlimited contacts on every plan including free, 300 emails/day free (~9,000/month), double opt-in forms, unsubscribe, campaigns, transactional API and SMTP included; monthly paid plans from about $9 for the one month you mail a large list, then back to free. Confirm current numbers on their pricing page. It replaces the Worker + KV + Turnstile + Resend build entirely. *Not Kit, MailerLite, Loops, Beehiiv, Buttondown or Mailchimp — all priced per contact, so a list you rarely email costs money every month. Not SendGrid (free tier gone, $19.95/mo). Amazon SES ($0.10 per 1,000, no contact pricing) with the list in Cloudflare KV is the fallback if Brevo ever changes terms: cheapest at any scale, most to run.*
- **[YOU]** Do not gate the sample report behind an email address: link the PDF directly on the site. The list is then purely "notify me of updates" — fewer people join, nothing has to be sent promptly to anyone, and the privacy obligations shrink to a newsletter's.
- **[YOU]** The sending address: `hello@openbiota.com`, verified in the email service (SPF, DKIM, DMARC records at Cloudflare), with replies forwarded to your existing mailbox by Email Routing. No Google Workspace unless you want a separate real inbox ($7/user/mo). Not `no-reply@`: a reply-able address is better for trust and deliverability.
- **[YOU]** The public sample report: your own report (you can consent to your own; the id is redacted to `SAMPLE2_A02` already) or a synthetic one from `openbiota simulate`. Real reads better; synthetic avoids the question entirely.
- **[YOU]** Copyright holder's legal name in `LICENSE` (currently "Erick Miller (OpenBiota)"). An entity later means re-licensing with contributors in the tree.
- **[YOU]** Paid product at launch, or free + interest list only? Section 9 is dead weight until this is answered.
- **[ME, done]** Repository `openbiota-org/openbiota-gut-health-test`, private, single squashed release commit authored as Erick Miller, no personal data in tree or history; `.gitignore` covers reads, results, references, third-party reports, secrets; no fixture derived from a real sample is tracked. — 2026-09-30.
- **[ME, done]** Report self-agreement: one judge of levels (`organisms.verdict`), one account of what is present (pathogen screen reconciled with the inventory, biofilm and disease tables name their catalogue and print the pooled reading), levels and deviations never point opposite ways; all held by tests that run on any sample. — 2026-09-30.

## 1. Make the repository public (this week; nothing else depends on the site)

- **[ME]** README status block at the top: measured recall/precision against the spec targets, the modules still labelled research-grade (mycobiome support rule, Myco-Score, determinant panels, pathogen rule, microbiome age, disease patterns), hardware and disk needed, and "research use only, not a diagnostic test". A visitor reads the limits before the claims.
- **[ME]** Pre-flight: `git grep` for secrets, tokens, absolute paths and local ids (done once; repeat at the moment of flipping); confirm `LICENSE`, `CITATION.cff` (new), `CODE_OF_CONDUCT.md`, `CONTRIBUTING.md`, `SECURITY.md`, issue and PR templates are in place.
- **[YOU]** Flip to public: `gh repo edit openbiota-org/openbiota-gut-health-test --visibility public --accept-visibility-change-consequences`. Then in Settings: branch protection on `main` requiring the `tests` check; enable Dependabot alerts + security updates, secret scanning + push protection (all free on public repos); enable Discussions as the support channel; set description, topics, website; upload the social preview image (`web/openbiota.com/img/openbiota-ogimg.png`). Delete the stray `erickmiller/gut-health-metagenomic-screen`.
- **[ME]** `CODEOWNERS`, Dependabot config (pip + GitHub Actions, weekly), and a `CHANGELOG.md` starting at v0.8.4.
- **[ME]** Tag `v0.8.4` with release notes; attach the small derived artefacts (section 4). Connect GitHub to Zenodo so every release gets a DOI automatically; cite it from README and the report's technical page.

## 2. Legal and policy (before the first email address is collected)

- **[ME]** Draft the privacy policy and terms pages (static HTML in `web/openbiota.com/`): what the site collects (an email address, via the named email service as processor; cookieless analytics), why, retention, how to unsubscribe and be deleted, contact `privacy@openbiota.com`, research-use-only statement, no medical advice. **[YOU]** review before publishing.
- **[ME]** Footer links to privacy, terms, contact and the research-use statement on every page.
- **[YOU]** DCO (simplest) or CLA before the first outside contribution. PolyForm Noncommercial with reserved commercial rights only works if contributors grant those rights; a DCO sign-off line in `CONTRIBUTING.md` is enough to start.
- **[YOU]** Decide the research-use-only posture in writing; the report says it everywhere, the site must too.
- **[ME, done]** Homepage wording: "free and open — source available under PolyForm Noncommercial", not "open source". — 2026-09-29.
- **[ME, done]** Third-party attribution (`docs/ATTRIBUTION.md`, generated from the reference locks and tool pins). — 2026-09-29.

## 3. Site and interest list (Cloudflare + the email service)

- **[YOU]** Cloudflare account; add `openbiota.com`; point nameservers (or transfer); enable Email Routing with `hello@`, `privacy@`, `security@` forwarding to your mailbox; enable Web Analytics (cookieless).
- **[ME]** Replace the GitHub Pages deploy with Cloudflare Pages: connect the repository in the Cloudflare dashboard (build command `mkdocs build --strict && python scripts/build_site.py --no-docs`, output `web/openbiota.com`), or swap the deploy step for `wrangler-action`. Remove the `DEPLOY_PAGES` gate and the Pages job; keep the strict build as the pull-request check. Custom domain + HTTPS on Pages. Add the analytics snippet.
- **[YOU]** Brevo account; verify `openbiota.com` as the sending domain (add the SPF/DKIM/DMARC records it gives you at Cloudflare); create one "updates" list with double opt-in; write nothing else — Brevo handles confirmation, unsubscribe and the suppression list.
- **[ME]** Point the site's existing dialog at Brevo's form endpoint (keep our markup, the honeypot and the live region; drop `emailEndpoint`/Worker assumptions from `script.js`); the confirmation email and the first "what's new" template, written in Brevo's editor.
- **[ME]** The public sample report once decided (section 0): redacted, read end to end against a checklist, linked directly from the site as a PDF — no signup in front of it.
- **[ME]** Verify every claim on the page against current output (counts, catalogue sizes, page counts, section numbers) and tone the hero copy to the measured capability. Add the privacy/terms links. Accessibility and performance pass (alt text, contrast, keyboard path through the dialog, Lighthouse). Test on mobile.
- **[ME, done]** OG/Twitter cards, favicon set, `robots.txt`, root `sitemap.xml`, 404 page — `scripts/build_site.py`. — 2026-09-29.

## 4. Reference data: an install anyone can reproduce

`refs/` is **1.3 TB** on the development machine and is three different things:
~960 GB of **upstream** databases that are not ours to ship (GlobDB 545, mycobiome
211, MetaPhlAn 85, sylph 48, SingleM 20, mOTUs 14, UHGG Kraken 16, host/decoys/
genomes ~18); ~145 GB of **cohort reads** downloaded once to build 14 MB of
cohort files that are already in git; and ~135 GB **derived by us** (pathogen
bundle 120 — competitive index 82, sequences 10, k-mer index 9, decoy k-mers 9,
assets 10 — rescue-panel Kraken DB 13, DIAMOND DBs 5) plus ~100 MB of small
derived files. We host none of the upstream data and never use S3 (egress on a
135 GB package is ~$12 per download).

- **[ME]** Tier A, in git (~20 MB, done): curated YAML, age models, `reference_ranges.json`, cohort manifests and files, pathogen catalogue, resolution cache, crosswalk summaries.
- **[ME]** Tier B lockfile `refs/MANIFEST.lock`: every upstream artefact with source URL, accession, byte size, SHA-256, licence, retrieval date — including the sylph GTDB database, MetaPhlAn Jan26 and the mycobiome set (`refs/mycobiome/reference_lock.json` already hashes those).
- **[ME]** `openbiota fetch-refs`: resumable, digest-verified download of Tier B from the publishers; `--verify-only`; `--build` to regenerate the derived artefacts from the existing scripts (pathogen bundle, rescue panel, DIAMOND DB) — the default path, zero hosting; `--download` to take the Zenodo snapshot instead where one exists.
- **[ME]** Install profiles with measured disk: *standard* (no GlobDB genome archive — confirmation fetches genomes on demand into `refs/genomes`) and *full*. Fix `docs/INSTALL.md`, which says ~300 GB.
- **[YOU]** Zenodo account (GitHub login) and a dataset record for the derived snapshot that fits a 50 GB record: k-mer index, decoy k-mers, rescue panel, DIAMOND DBs, small derived files (~36 GB). The 82 GB competitive index stays rebuild-only unless people ask; then Cloudflare R2 (~$1.30/month at 85 GB, zero egress) on the account you already have.
- **[ME]** Cite the dataset DOI from README, `docs/INSTALL.md` and the report's technical page; a scheduled CI job re-verifies every lockfile URL and digest so a moved upstream file is caught before a user hits it.
- **[ME]** `openbiota doctor` tells a new user exactly what is missing and the command that fixes it, lane by lane.
- **[ME]** `make demo`: `openbiota simulate` on a small synthetic community against the standard profile, end to end to a PDF, so a visitor sees the pipeline work before committing a terabyte.

## 5. Validation still owed before the numbers are claims

- **[ME]** Share accuracy on mock communities: the lane benchmark measures detection (94.2% recall, 92.6% precision); nothing yet measures how far a split, estimated or absorbed share lands from a known abundance. Add the abundance-error measure to `scripts/benchmark_lanes.py` and print it in `docs/MEASURED_CAPABILITY.md`.
- **[ME]** Widen the lane cohort: the expansion-only organisms are ranked against 100 prefix-sampled adults, which is coarse. More ENA samples cost only compute.
- **[ME]** Mycobiome: qualify the support rule and discovery threshold on the Avershina 2025 mock workflow plus simulated 0.0001–1% fungal fractions (the "starting values, not yet qualified on mocks" label stays until then); run FungiGut unmodified as the reproduction baseline; strain validation on downsampled isolate reads and mixtures; the determinant panel (ERG11/FKS, ECE1) with species-specific coordinates; donor–recipient strain comparison; rebuild the species-discovery index split so it runs in ~30 GB.
- **[YOU]** Whether the Myco-Score (weights are policy, labelled Experimental) appears on the public sample report before the mock qualification lands.
- **[ME]** Modulator-registry coverage test: a flagged organism in any run with no registry entry fails, so a new organism forces a curated entry. **[YOU]** review cadence for its citations (suggest: each reference release).

## 6. Software release hygiene

- **[ME, done]** CI runs ruff and the suite on Python 3.11–3.13 on every push; data-dependent tests skip on a clean checkout. — 2026-09-29.
- **[ME]** Pin and lock Python dependencies (`requirements.lock` covers only the strain tools today).
- **[ME]** Report versioning: `make bump-patch|minor|major`, `CHANGELOG.md`, the report-surface test that fails when the report gains or loses a reading without a version bump, and `reference_release` as a version string with the rule for when a reference change forces a bump.
- **[ME]** Release workflow on tag: run the suite, build the small derived artefacts, attach them, publish notes; Zenodo DOI via the integration.
- **[YOU]** PyPI or not. If yes, register `openbiota` and store a token as a repository secret.

## 7. Documentation

- **[ME, done]** `docs/` current for the merged inventory, detection lanes, one share per organism, level against the reference, attention rule, pathogen reconciliation. — 2026-09-30.
- **[ME]** "Interpreting your report" guide — the most-requested thing a reader will want and the target of the second email.
- **[ME]** Document `--strict-cdiff` and every other non-default interpretation rule in `docs/USAGE.md`.
- **[ME, done]** Hardware and runtime expectations measured from the six runs — `docs/INSTALL.md` (disk figure to be corrected under section 4). — 2026-09-29.

## 8. Before the first real reader

- **[ME]** Read one complete report end to end, all ~330 pages, against a checklist.
- **[YOU]** Have a clinician read the report and confirm nothing reads as a diagnosis.
- **[YOU]** Decide what happens when a report finds something genuinely alarming. There is no protocol today, and the pathogen section can produce one.

## 9. If selling a product (skip entirely if free-only at launch)

- **[YOU]** Entity, business bank account, accountant conversation about sales-tax nexus.
- **[YOU]** Stripe account; price; whether the kit is bundled.
- **[YOU]** Sequencing partner, per-sample pricing and turnaround — the long pole.
- **[YOU]** Collection kits: supplier, labelling, UN3373 shipping category, return postage.
- **[YOU]** CLIA/CAP posture in writing. Research-use-only avoids it; marketing a report as informing a health decision changes the analysis.
- **[YOU]** Liability insurance and a support channel with a real SLA.
- **[ME]** Order flow (Stripe Checkout, order record, kit-dispatch trigger, tracking email); sample intake (barcode → order, chain of custody, status page); secure results delivery (authenticated download, expiring links, encryption at rest, deletion policy, data export); pipeline automation (queue, retries, alerting, per-sample cost ledger); capacity plan (~35 min/sample on one box today).

---

### Cost summary

| | Account | Monthly cost | Ceiling before anything costs money |
|---|---|---|---|
| Code, CI, releases, issues, discussions | GitHub (have) | $0 | public repo: unlimited |
| DNS, hosting, inbound mail, analytics, future backend | Cloudflare | $0 | Pages unlimited bandwidth; Workers 100k requests/day; R2 10 GB |
| Interest list, double opt-in, campaigns, sending domain | Brevo | $0 | unlimited contacts; 300 emails/day — pay (~$9+) only in a month you mail more than that |
| Citable data and software snapshots | Zenodo | $0 | 50 GB per record |
| Domain | registrar | ~$10/yr | — |
