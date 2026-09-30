/* OpenBiota — the research catalog page.
 * Reads research-data.js (generated from the project research catalog) and
 * renders a searchable table. No build step, framework, or server required;
 * research.html can be opened straight from disk.
 */
(() => {
'use strict';

const data = window.OPENBIOTA_RESEARCH;
const body = document.querySelector('[data-research-body]');
const statusText = document.querySelector('[data-research-status]');
if (!body) return;
if (!data || !Array.isArray(data.records)) {
  statusText.textContent = 'The catalog data could not be loaded. Keep research-data.js beside research.html.';
  return;
}

// Field positions in each record array, as written by the generator.
const ID = 0, NAME = 1, YEAR = 2, WHAT = 3, HOW = 4, INST = 5, COUNTRY = 6,
      BASIS = 7, TYPE = 8, FAMILY = 9, LINK = 10, EXTRA = 11, NOTES = 12,
      SPECS = 13, DOCS = 14, ACCESSION = 15, RESOURCE = 16;

const PAGE_SIZE = 30, PAGE_STEP = 60;
const BASIS_TEXT = {
  S: 'Documented study, sample, or experiment location.',
  A: 'Location of the principal author institutions; the study setting is not independently established.',
  M: 'Location of the organization that develops or issues the resource.'
};
const HOSTS = {
  'doi.org': 'DOI', 'europepmc.org': 'Europe PMC', 'pubmed.ncbi.nlm.nih.gov': 'PubMed',
  'pmc.ncbi.nlm.nih.gov': 'PubMed Central', 'ncbi.nlm.nih.gov': 'NCBI',
  'ftp.ncbi.nlm.nih.gov': 'NCBI', 'ebi.ac.uk': 'EMBL-EBI', 'ftp.ebi.ac.uk': 'EMBL-EBI',
  'github.com': 'GitHub', 'bitbucket.org': 'Bitbucket', 'zenodo.org': 'Zenodo',
  'cdc.gov': 'CDC', 'wwwnc.cdc.gov': 'CDC', 'fda.gov': 'FDA', 'nih.gov': 'NIH',
  'ods.od.nih.gov': 'NIH', 'niddk.nih.gov': 'NIH', 'hhs.gov': 'HHS',
  'nature.com': 'Nature', 'science.org': 'Science', 'cell.com': 'Cell',
  'academic.oup.com': 'Oxford Academic', 'gut.bmj.com': 'BMJ Gut', 'bmj.com': 'BMJ',
  'sciencedirect.com': 'ScienceDirect', 'frontiersin.org': 'Frontiers',
  'journals.plos.org': 'PLOS', 'elifesciences.org': 'eLife', 'mdpi.com': 'MDPI',
  'biorxiv.org': 'bioRxiv', 'medrxiv.org': 'medRxiv', 'clinicaltrials.gov': 'ClinicalTrials.gov',
  'ega-archive.org': 'EGA', 'springer.com': 'Springer', 'link.springer.com': 'Springer',
  'media.springernature.com': 'Springer Nature', 'ecfr.gov': 'eCFR',
  'who.int': 'WHO', 'wiley.com': 'Wiley', 'onlinelibrary.wiley.com': 'Wiley',
  'tandfonline.com': 'Taylor & Francis', 'journals.asm.org': 'ASM', 'jci.org': 'JCI',
  'pnas.org': 'PNAS', 'oup.com': 'Oxford Academic', 'figshare.com': 'figshare'
};

const numberFormat = new Intl.NumberFormat('en-US');
const number = value => numberFormat.format(value);
const escapeHtml = value => String(value).replace(/[&<>"']/g,
  character => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[character]));
// Notes often cite an erratum or the original article; make those reachable.
const linkifyNote = text => escapeHtml(text).replace(/https?:\/\/[^\s<>"')\]]+/g, match => {
  const trailing = /[.,;:]+$/.exec(match);
  const url = trailing ? match.slice(0, -trailing[0].length) : match;
  return `<a href="${url}" target="_blank" rel="noopener noreferrer">${url}</a>${trailing ? trailing[0] : ''}`;
});

function expand(code) {
  if (!code) return '';
  const prefix = data.prefixes[code[0]];
  return prefix ? prefix + code.slice(1) : code.slice(1);
}
// Titles carry typographic primes and dashes (2′-FL, Crohn’s, en-dashes).
// Folding them on both sides keeps an ASCII search from silently missing.
const fold = value => value.toLowerCase()
  .replace(/[\u2018\u2019\u02bc\u2032\u0060']/g, '')
  .replace(/[\u2013\u2014\u2212]/g, '-');

function hostName(url) {
  const match = /^https?:\/\/([^/]+)/.exec(url);
  if (!match) return 'Source';
  const host = match[1].replace(/^www\./, '');
  return HOSTS[host] || HOSTS[host.split('.').slice(-2).join('.')] || host;
}

// One lowercase haystack per record keeps typing responsive across 2,000 rows.
const rows = data.records.map(record => {
  const institutions = record[INST].map(index => data.institutions[index]);
  const countries = record[COUNTRY].map(index => data.countries[index]);
  // The catalogs describe study design in their own words; `family` is the
  // coarse grouping the filter offers, `type` is the exact source wording.
  const type = data.types[record[TYPE]];
  const family = data.families[record[FAMILY]];
  return {
    record, institutions, countries, type, family,
    link: expand(record[LINK]),
    extra: record[EXTRA].map(expand),
    haystack: fold([record[ID], record[NAME], record[YEAR], record[WHAT], record[HOW],
                    record[ACCESSION], type, family, institutions.join(' '),
                    countries.join(' '), record[NOTES].join(' ')].join(' '))
  };
});

// ------------------------------------------------------------------ controls
const form = document.querySelector('.research-controls');
const searchInput = document.querySelector('#research-search-input');
const clearButton = document.querySelector('[data-research-clear]');
const kindSelect = document.querySelector('#research-kind');
const typeSelect = document.querySelector('#research-type');
const countrySelect = document.querySelector('#research-country');
const institutionSelect = document.querySelector('#research-institution');
const sortSelect = document.querySelector('#research-sort');
const moreButton = document.querySelector('[data-research-more]');
const progressText = document.querySelector('[data-research-progress]');
const emptyPanel = document.querySelector('[data-research-empty]');
const resetButtons = document.querySelectorAll('[data-research-reset]');

for (const [key, value] of Object.entries(data.counts || {})) {
  for (const target of document.querySelectorAll(`[data-research-count="${key}"]`)) {
    target.textContent = number(value);
  }
}

function fillSelect(select, entries) {
  const fragment = document.createDocumentFragment();
  for (const [value, label, count] of entries) {
    const option = document.createElement('option');
    option.value = value;
    option.textContent = count ? `${label} (${number(count)})` : label;
    fragment.append(option);
  }
  select.append(fragment);
}
function tally(getValues) {
  const counts = new Map();
  for (const row of rows) {
    for (const value of getValues(row)) counts.set(value, (counts.get(value) || 0) + 1);
  }
  return [...counts.entries()]
    .sort((a, b) => a[0].localeCompare(b[0], 'en'))
    .map(([value, count]) => [value, value, count]);
}
fillSelect(typeSelect, tally(row => [row.family]));
fillSelect(countrySelect, tally(row => row.countries));
fillSelect(institutionSelect, tally(row => row.institutions));

// -------------------------------------------------------------------- filters
const state = { query: '', kind: 'all', type: 'all', country: 'all', institution: 'all', sort: 'grouped', shown: PAGE_SIZE };
const open = new Set();
let visible = rows;

function matchesKind(row) {
  const record = row.record;
  switch (state.kind) {
    case 'publications': return !record[RESOURCE];
    case 'resources': return Boolean(record[RESOURCE]);
    case 'trials': return ['Randomized trial', 'Review or meta-analysis',
                           'Guideline or clinical reference', 'Case report'].includes(row.family);
    case 'notes': return record[NOTES].length > 0;
    default: return true;
  }
}
function compare(a, b) {
  const first = a.record, second = b.record;
  const byName = first[NAME].localeCompare(second[NAME], 'en');
  switch (state.sort) {
    case 'name': return byName;
    case 'newest': return (Number(second[YEAR]) || 0) - (Number(first[YEAR]) || 0) || byName;
    case 'oldest': return (Number(first[YEAR]) || 9999) - (Number(second[YEAR]) || 9999) || byName;
    case 'id': return first[ID].localeCompare(second[ID], 'en');
    // Published research leads, then datasets, databases, and software.
    default: return first[RESOURCE] - second[RESOURCE] || byName;
  }
}
function applyFilters() {
  const terms = fold(state.query).split(/\s+/).filter(Boolean);
  visible = rows.filter(row => {
    if (!matchesKind(row)) return false;
    if (state.type !== 'all' && row.family !== state.type) return false;
    if (state.country !== 'all' && !row.countries.includes(state.country)) return false;
    if (state.institution !== 'all' && !row.institutions.includes(state.institution)) return false;
    return terms.every(term => row.haystack.includes(term));
  });
  visible.sort(compare);
}

// ------------------------------------------------------------------ rendering
function list(values, limit) {
  if (!values.length) return '';
  const shown = values.slice(0, limit).map(escapeHtml).join('<span class="research-sep">·</span>');
  const rest = values.length - limit;
  return rest > 0 ? `${shown}<span class="research-rest">+${rest} more</span>` : shown;
}

function detailHtml(row) {
  const record = row.record;
  const parts = [];
  const meta = [`<span><b>Catalog record</b>${escapeHtml(record[ID])}</span>`,
                `<span><b>Evidence type</b>${escapeHtml(row.type)}</span>`];
  if (record[YEAR]) meta.push(`<span><b>Year</b>${record[YEAR]}</span>`);
  if (record[ACCESSION]) meta.push(`<span><b>Accession</b>${escapeHtml(record[ACCESSION])}</span>`);
  if (record[SPECS].length) meta.push(`<span><b>Build specs</b>${record[SPECS].map(escapeHtml).join(', ')}</span>`);
  parts.push(`<div class="research-detail-meta">${meta.join('')}</div>`);

  if (record[NOTES].length) {
    parts.push(`<div class="research-detail-block"><h4>Provenance notes &amp; interpretation limits</h4><ul>${
      record[NOTES].map(note => `<li>${linkifyNote(note)}</li>`).join('')}</ul></div>`);
  }
  if (row.institutions.length > 2 || row.countries.length > 2) {
    const blocks = [];
    if (row.institutions.length > 2) blocks.push(`<p><b>All institutions:</b> ${row.institutions.map(escapeHtml).join('; ')}</p>`);
    if (row.countries.length > 2) blocks.push(`<p><b>All countries:</b> ${row.countries.map(escapeHtml).join('; ')}</p>`);
    parts.push(`<div class="research-detail-block">${blocks.join('')}</div>`);
  }
  if (BASIS_TEXT[record[BASIS]] && row.countries.length) {
    parts.push(`<p class="research-detail-basis">${escapeHtml(BASIS_TEXT[record[BASIS]])}</p>`);
  }
  if (record[DOCS].length) {
    parts.push(`<p class="research-detail-basis">Recorded in: ${
      record[DOCS].map(index => escapeHtml(data.documents[index])).join(', ')}.</p>`);
  }
  if (row.extra.length) {
    parts.push(`<div class="research-detail-links">${row.extra.map(url =>
      `<a href="${escapeHtml(url)}" target="_blank" rel="noopener noreferrer">${escapeHtml(hostName(url))}<svg class="icon icon-external" width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" focusable="false"><path d="M15 3h6v6M10 14 21 3M21 14v5a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V5a2 2 0 0 1 2-2h5"/></svg></a>`).join('')}</div>`);
  }
  return parts.join('');
}

function rowHtml(row) {
  const record = row.record;
  const id = record[ID];
  const expanded = open.has(id);
  const meta = [escapeHtml(row.family)];
  if (record[YEAR]) meta.push(record[YEAR]);
  if (record[ACCESSION]) meta.push(escapeHtml(record[ACCESSION]));
  const linkCell = row.link
    ? `<a class="research-link" href="${escapeHtml(row.link)}" target="_blank" rel="noopener noreferrer">${
        record[RESOURCE] ? 'Open' : 'Read'}<svg class="icon icon-external" width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" focusable="false"><path d="M15 3h6v6M10 14 21 3M21 14v5a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V5a2 2 0 0 1 2-2h5"/></svg></a><span class="research-link-host">${
        escapeHtml(hostName(row.link))}</span>`
    : '';
  const noteFlag = record[NOTES].length
    ? `<span class="research-flag" title="${record[NOTES].length} provenance note${record[NOTES].length > 1 ? 's' : ''}">note${record[NOTES].length > 1 ? 's' : ''}</span>`
    : '';

  return `<tr class="research-row${expanded ? ' is-open' : ''}" id="${escapeHtml(id)}">
<td class="col-name" data-label="Research"><div class="research-name">${escapeHtml(record[NAME])}</div><div class="research-meta">${
  meta.join('<span class="research-sep">·</span>')}${noteFlag}</div><button type="button" class="research-info" data-research-toggle="${
  escapeHtml(id)}" aria-expanded="${expanded}" aria-controls="detail-${escapeHtml(id)}"><span class="research-info-mark" aria-hidden="true">i</span><span>${
  expanded ? 'Hide details' : 'Details'}</span></button></td>
<td class="col-what" data-label="What it is">${escapeHtml(record[WHAT])}</td>
<td class="col-how" data-label="How it informed OpenBiota">${escapeHtml(record[HOW])}</td>
<td class="col-inst" data-label="Institutions">${list(row.institutions, 2)}</td>
<td class="col-country" data-label="Country">${list(row.countries, 2)}</td>
<td class="col-link" data-label="Source">${linkCell}</td>
</tr>
<tr class="research-detail" id="detail-${escapeHtml(id)}"${expanded ? '' : ' hidden'}><td colspan="6">${
  expanded ? detailHtml(row) : ''}</td></tr>`;
}

function render() {
  const slice = visible.slice(0, state.shown);
  body.innerHTML = slice.map(rowHtml).join('');
  const total = visible.length;
  const filtered = total !== rows.length || state.query;

  emptyPanel.hidden = total > 0;
  document.querySelector('.research-table-wrap').hidden = total === 0;
  moreButton.hidden = slice.length >= total;
  progressText.textContent = total > slice.length
    ? `Showing ${number(slice.length)} of ${number(total)} sources`
    : total === 1 ? '1 matching source shown'
    : total ? `All ${number(total)} matching sources shown` : '';

  statusText.textContent = filtered
    ? `${number(total)} matching source${total === 1 ? '' : 's'}`
    : `${number(rows.length)} sources — ${number(data.counts.publications)} published studies and ${number(data.counts.resources)} datasets, databases, and software resources`;
  for (const button of resetButtons) button.hidden = !filtered;
}

function update({ resetPage = true } = {}) {
  if (resetPage) state.shown = PAGE_SIZE;
  applyFilters();
  render();
}

// --------------------------------------------------------------- interaction
let debounce;
searchInput.addEventListener('input', () => {
  clearButton.hidden = !searchInput.value;
  clearTimeout(debounce);
  debounce = setTimeout(() => {
    state.query = searchInput.value.trim();
    update();
  }, 120);
});
searchInput.addEventListener('search', () => {
  state.query = searchInput.value.trim();
  clearButton.hidden = !searchInput.value;
  update();
});
clearButton.addEventListener('click', () => {
  searchInput.value = '';
  state.query = '';
  clearButton.hidden = true;
  searchInput.focus();
  update();
});
form.addEventListener('submit', event => {
  event.preventDefault();
  state.query = searchInput.value.trim();
  update();
});
for (const [select, key] of [[kindSelect, 'kind'], [typeSelect, 'type'],
                             [countrySelect, 'country'], [institutionSelect, 'institution'],
                             [sortSelect, 'sort']]) {
  select.addEventListener('change', () => {
    state[key] = select.value;
    update();
  });
}
moreButton.addEventListener('click', () => {
  const firstNew = state.shown;
  state.shown += PAGE_STEP;
  render();
  body.children[firstNew * 2]?.querySelector('.research-info')?.focus({ preventScroll: true });
});
for (const button of resetButtons) {
  button.addEventListener('click', () => {
    searchInput.value = '';
    clearButton.hidden = true;
    Object.assign(state, { query: '', kind: 'all', type: 'all', country: 'all', institution: 'all' });
    kindSelect.value = typeSelect.value = countrySelect.value = institutionSelect.value = 'all';
    update();
    searchInput.focus();
  });
}
body.addEventListener('click', event => {
  const toggle = event.target.closest('[data-research-toggle]');
  if (!toggle) return;
  const id = toggle.dataset.researchToggle;
  const row = visible.find(entry => entry.record[ID] === id);
  const detail = document.querySelector(`#detail-${CSS.escape(id)}`);
  if (!row || !detail) return;
  const expanded = !open.has(id);
  if (expanded) open.add(id); else open.delete(id);
  toggle.setAttribute('aria-expanded', String(expanded));
  toggle.querySelector('span:last-child').textContent = expanded ? 'Hide details' : 'Details';
  toggle.closest('tr').classList.toggle('is-open', expanded);
  detail.firstElementChild.innerHTML = expanded ? detailHtml(row) : '';
  detail.hidden = !expanded;
});

// A catalog ID in the address (research.html#P0031) opens that record.
function openFromHash() {
  // P0031 from the historical catalog, P0830009 / R083… / V083… / I083… from v0.8.3.
  const id = decodeURIComponent(location.hash.slice(1)).toUpperCase();
  if (!/^[PRVI]\d{4,7}$/.test(id) || !rows.some(row => row.record[ID] === id)) return;
  searchInput.value = id;
  clearButton.hidden = false;
  state.query = id;
  open.add(id);
  update();
  document.getElementById(id)?.scrollIntoView({ block: 'center' });
}

const initialQuery = new URLSearchParams(location.search).get('q');
if (initialQuery) {
  searchInput.value = initialQuery;
  state.query = initialQuery.trim();
  clearButton.hidden = false;
}
update();
openFromHash();
window.addEventListener('hashchange', openFromHash);

})();
