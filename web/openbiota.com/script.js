/* OpenBiota — portable, framework-free landing page.
 * Open index.html directly or upload this folder to any static web host.
 * No build step, npm, platform account, or OpenAI server is required.
 *
 * Optional setup: replace repositoryUrl and emailEndpoint below.
 * emailEndpoint accepts a JSON POST and returns a JSON response.
 * Keep mail-provider API keys on your server, never in this file.
 * See README.md for the small request/response contract.
 */
const OPENBIOTA = {
  repositoryUrl: 'https://github.com/openbiota/OpenBiota-Gut-Health-Report',
  emailEndpoint: '',
  eventsEndpoint: ''
};

(() => {
'use strict';

// Navigation and expandable sections.
const toggle=document.querySelector('.menu-toggle');
const navigation=document.querySelector('#main-nav');
function closeNavigation(){toggle?.setAttribute('aria-expanded','false');toggle?.setAttribute('aria-label','Open navigation');navigation?.classList.remove('open');}
toggle?.addEventListener('click',event=>{const open=toggle.getAttribute('aria-expanded')!=='true';toggle.setAttribute('aria-expanded',String(open));toggle.setAttribute('aria-label',open?'Close navigation':'Open navigation');navigation.classList.toggle('open',open);if(open&&event.detail===0)navigation.querySelector('a')?.focus();});
document.addEventListener('keydown',event=>{if(event.key==='Escape'&&navigation?.classList.contains('open')){closeNavigation();toggle.focus();}});
// Repeated activation must also reveal a destination that was collapsed again.
function revealAnchor(hash=location.hash){
  let target;try{target=document.getElementById(decodeURIComponent(hash.slice(1)));}catch{return;}
  if(!target)return;let parent=target.closest('details');while(parent){parent.open=true;parent=parent.parentElement.closest('details');}
  requestAnimationFrame(()=>target.scrollIntoView({block:'start'}));
}
document.addEventListener('click',event=>{
  const anchor=event.target.closest('a');
  if(!event.target.closest('.header')||anchor?.closest('.header'))closeNavigation();
  if(event.defaultPrevented)return;
  if(anchor?.getAttribute('href')?.startsWith('#'))revealAnchor(anchor.getAttribute('href'));
});
window.addEventListener('hashchange',()=>revealAnchor());
if(location.hash&&!['#get-report','#sample-report'].includes(location.hash))revealAnchor();

// The wordmark returns to the top of the home page and leaves the address bar
// clean — no "#top" — as if the page had just been opened. Modified clicks
// (new tab, new window) keep their default behaviour.
const brand=document.querySelector('.header .brand');
brand?.addEventListener('click',event=>{
  if(!brand.getAttribute('href').startsWith('#')||event.button!==0||event.metaKey||event.ctrlKey||event.shiftKey||event.altKey)return;
  if(location.hash){try{history.pushState(history.state,'',location.pathname+location.search);}catch{return;/* file:// may refuse history changes; "#top" still works. */}}
  event.preventDefault();window.scrollTo({top:0});
});

// Email modal and optional demand measurement.
const knownSources=new Set(['direct','hackernews','producthunt','x','linkedin','instagram','tiktok','youtube','facebook','reddit','newsletter','influencer','other']);
const sourceAliases={hn:'hackernews',ph:'producthunt',twitter:'x'};
const cleanCampaign=value=>String(value||'').toLowerCase().replace(/[^a-z0-9_-]/g,'').slice(0,64);
function currentSource(){
  const params=new URLSearchParams(location.search),tag=params.get('utm_source')?.toLowerCase();
  if(tag){const source=sourceAliases[tag]||tag;return {source:knownSources.has(source)?source:'other',campaign:cleanCampaign(params.get('utm_campaign'))};}
  try{const ref=new URL(document.referrer);if(ref.origin===location.origin)return null;const host=ref.hostname.replace(/^www\./,'');const domains={'news.ycombinator.com':'hackernews','producthunt.com':'producthunt','x.com':'x','t.co':'x','twitter.com':'x','linkedin.com':'linkedin','instagram.com':'instagram','tiktok.com':'tiktok','youtube.com':'youtube','facebook.com':'facebook','reddit.com':'reddit'};const domain=Object.keys(domains).find(d=>host===d||host.endsWith('.'+d));return {source:domain?domains[domain]:'other',campaign:''};}catch{return null;}
}
let attribution=currentSource();
try{if(!attribution){const stored=JSON.parse(sessionStorage.getItem('openbiota-attribution')||'null');if(stored&&knownSources.has(stored.source))attribution={source:stored.source,campaign:cleanCampaign(stored.campaign)};}attribution??={source:'direct',campaign:''};sessionStorage.setItem('openbiota-attribution',JSON.stringify(attribution));}catch{attribution??={source:'direct',campaign:''};}
// Optional analytics: disabled until you supply an endpoint.
// Also emits a local browser event so you can attach your preferred analytics.
function track(event) {
  const detail = { event, ...attribution };
  window.dispatchEvent(new CustomEvent('openbiota:interest', { detail }));
  if (!OPENBIOTA.eventsEndpoint) return;
  fetch(OPENBIOTA.eventsEndpoint, {
    method: 'POST',
    headers: { 'Content-Type': 'application/json' },
    body: JSON.stringify(detail),
    keepalive: true
  }).catch(() => {});
}
track('landing_view');

const reportDialog=document.querySelector('#report-dialog');
const interestForm=document.querySelector('#report-interest-form');
const formView=reportDialog.querySelector('.request-form-view');
const successView=reportDialog.querySelector('.request-success');
const statusText=document.querySelector('#interest-status');
const emailInput=document.querySelector('#report-email');
const submitButton=interestForm.querySelector('button[type="submit"]');
let returnFocus=null,savedScroll=0,bodyStyles=null,submitting=false;

// Use your own repository URL. No redirect service is required.
for (const link of document.querySelectorAll('[data-openbiota-code]')) {
  link.href = OPENBIOTA.repositoryUrl;
  link.target = '_blank';
  link.rel = 'noopener noreferrer';
}

function updateDeliveryNote() {
  document.querySelector('#email-note').textContent = OPENBIOTA.emailEndpoint
    ? 'The software is 100% free and open source.'
    : 'The software is 100% free and open source.';
}
function openReportRequest(kind='sample',placement='direct',trigger=null){
  if(reportDialog.open)return;
  closeNavigation();returnFocus=trigger||document.activeElement;savedScroll=window.scrollY;
  if(!submitting){
  interestForm.elements.interest_kind.value=kind==='own_report'?'own_report':'sample';
  interestForm.elements.placement.value=placement;
  document.querySelector('#request-title').textContent=kind==='own_report'?'Your microbiome report.  100% free, open source.':'See what the report reveals. 100% free, open source.';
  document.querySelector('#request-description').textContent=kind==='own_report'
    ?'Preview a free example report and get the guide to generating your own gut health report. OpenBiota is free, open-source software you can use, explore, and build on.'
    :'Preview a free example report and get the guide to generating your own. OpenBiota is free open-source software, shared freely for anyone to use and improve.';
  formView.hidden=false;successView.hidden=true;statusText.textContent='';statusText.classList.remove('error');emailInput.removeAttribute('aria-invalid');
  reportDialog.setAttribute('aria-labelledby','request-title');reportDialog.setAttribute('aria-describedby','request-description');
  updateDeliveryNote();
  }
  bodyStyles={position:document.body.style.position,top:document.body.style.top,width:document.body.style.width,overflow:document.body.style.overflow,paddingRight:document.body.style.paddingRight};
  const scrollbarWidth=window.innerWidth-document.documentElement.clientWidth;if(scrollbarWidth>0)document.body.style.paddingRight=`${scrollbarWidth}px`;
  document.body.style.position='fixed';document.body.style.top=`-${savedScroll}px`;document.body.style.width='100%';document.body.style.overflow='hidden';
  reportDialog.showModal();reportDialog.scrollTop=0;emailInput.focus({preventScroll:true});
  if(!submitting)track(kind==='own_report'?'own_report_open':'sample_open');
}
document.querySelectorAll('[data-report-request]').forEach(button=>button.addEventListener('click',()=>openReportRequest(button.dataset.reportRequest,button.dataset.placement,button)));
document.querySelectorAll('[data-close-dialog]').forEach(button=>button.addEventListener('click',()=>reportDialog.close()));
reportDialog.addEventListener('click',event=>{if(event.target!==reportDialog)return;const r=reportDialog.getBoundingClientRect();if(event.clientX<r.left||event.clientX>r.right||event.clientY<r.top||event.clientY>r.bottom)reportDialog.close();});
reportDialog.addEventListener('close',()=>{
  if(bodyStyles)Object.assign(document.body.style,bodyStyles);
  const oldBehavior=document.documentElement.style.scrollBehavior;document.documentElement.style.scrollBehavior='auto';window.scrollTo(0,savedScroll);returnFocus?.focus({preventScroll:true});document.documentElement.style.scrollBehavior=oldBehavior;
  if(['#get-report','#sample-report'].includes(location.hash)){try{const url=new URL(location.href);url.hash='';history.replaceState(null,'',url.href);}catch{/* Some browsers restrict file:// history changes. */}}
});
emailInput.addEventListener('input',()=>{emailInput.removeAttribute('aria-invalid');emailInput.setCustomValidity('');if(!submitting){statusText.textContent='';statusText.classList.remove('error');}});
emailInput.addEventListener('invalid',()=>emailInput.setAttribute('aria-invalid','true'));
interestForm.addEventListener('submit',async event=>{
  event.preventDefault();if(submitting)return;
  emailInput.value=emailInput.value.trim();
  if(!/^[^\s@\x00-\x1f]+@[^\s@\x00-\x1f]+\.[^\s@\x00-\x1f]+$/.test(emailInput.value))emailInput.setCustomValidity('Please enter a valid email address.');
  if(!interestForm.reportValidity())return;
  if(!OPENBIOTA.emailEndpoint){
    statusText.classList.add('error');
    statusText.textContent='Email requests aren’t available yet. Please check back soon.';
    return;
  }
  submitting=true;submitButton.disabled=true;interestForm.setAttribute('aria-busy','true');submitButton.querySelector('span').textContent='Sending your request…';statusText.textContent='';statusText.classList.remove('error');
  const controller=new AbortController(),timeout=setTimeout(()=>controller.abort(),20000);
  try{
    const response=await fetch(OPENBIOTA.emailEndpoint,{method:'POST',headers:{'Content-Type':'application/json',Accept:'application/json'},body:JSON.stringify({...Object.fromEntries(new FormData(interestForm)),...attribution}),signal:controller.signal});
    const data=await response.json();
    if(!response.ok||data.ok===false||data.error||data.errors)throw new Error(data.error||'Your request could not be saved. Please try again.');
    const sent=data.delivery==='sent';
    document.querySelector('#success-title').textContent=sent?'Check your inbox.':'Your request is received.';
    document.querySelector('#success-message').textContent=sent
      ?'We’ve emailed your free sample report and setup instructions. If you don’t see it shortly, check your spam folder.'
      :'Thank you—your request for the free sample report and setup instructions has been received.';
    track('report_request_submitted');
    formView.hidden=true;successView.hidden=false;reportDialog.setAttribute('aria-labelledby','success-title');reportDialog.setAttribute('aria-describedby','success-message');reportDialog.scrollTop=0;
    if(reportDialog.open)successView.focus({preventScroll:true});
    interestForm.reset();
  }catch(error){statusText.classList.add('error');statusText.textContent=error.name==='AbortError'?'The request took too long. Please try again in a moment.':error.message||'Your request could not be saved. Please try again.';}
  finally{clearTimeout(timeout);submitting=false;submitButton.disabled=false;interestForm.removeAttribute('aria-busy');submitButton.querySelector('span').textContent='Email me the free report';}
});
function openRequestHash(){if(location.hash==='#get-report')openReportRequest('own_report','shared-link');if(location.hash==='#sample-report')openReportRequest('sample','shared-link');}
window.addEventListener('hashchange',openRequestHash);openRequestHash();

// Section URLs. As the home page scrolls, the address bar names the section
// under the middle of the viewport, so any part of the page can be shared or
// reloaded; above the first section it is the bare home page. Writes wait
// for scrolling to settle and happen only on change: Safari refuses more
// than 100 history updates in 30 seconds, and replaceState never scrolls or
// adds history entries, so the back button still follows the reader's clicks.
const sectionRoot=document.querySelector('main[data-section-urls]');
const sections=sectionRoot?[...sectionRoot.querySelectorAll(':scope>section[id]')]:[];
function replaceHash(hash){
  if((location.hash||'')===hash)return;
  try{history.replaceState(history.state,'',location.pathname+location.search+hash);}catch{/* file:// may refuse history changes. */}
}
function syncSectionUrl(){
  if(reportDialog.open)return;// The open modal pins the body, so positions are meaningless until it closes.
  const middle=window.innerHeight/2;let current=null;
  for(const section of sections){if(section.getBoundingClientRect().top<=middle)current=section;else break;}
  if(!current){replaceHash('');return;}
  // A deep link into the section — an expandable, an article — is more precise than the section itself; keep it until the reader moves on.
  let target=null;try{target=document.getElementById(decodeURIComponent(location.hash.slice(1)));}catch{/* Malformed hash: replace it. */}
  if(target&&target!==current&&current.contains(target))return;
  replaceHash('#'+current.id);
}
if(sections.length){
  let settle=0;const schedule=()=>{clearTimeout(settle);settle=setTimeout(syncSectionUrl,150);};
  window.addEventListener('scroll',schedule,{passive:true});
  window.addEventListener('resize',schedule);
  schedule();
}

})();
