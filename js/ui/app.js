// BublaGene: application state and wiring.

import { h, clear, fmt, toast } from './dom.js';
import { circularMap, linearMap } from './map.js';
import { sequenceView } from './seqview.js';
import { showMenu, closeMenu } from './menu.js';
import { saveSession, loadSession, clearSession, timeAgo } from './session.js';
import { readAlignments } from '../formats/bam.js';
import { readTraces } from '../formats/ztr.js';
import { traceView } from './traceview.js';
import { PANELS, renderPanel } from './panels.js';
import { exampleDocument } from './sample.js';
import { readDna, writeDna } from '../formats/dna.js';
import { readGenBank, writeGenBank, readFasta, writeFasta } from '../formats/genbank.js';
import { newDocument, cloneDoc, replaceRange, rotateOrigin, reverseComplementDoc, addHistory, featureSpan, featureSequence, codingSequence, newId, colorFor } from '../model.js';
import { findSites, ENZYME_BY_NAME, COMMON_SET } from '../enzymes.js';
import { findBindingSites } from '../primers.js';
import { gcContent, meltingTemp, revcomp, subseq, rangeLength, cleanSequence, translate } from '../seq.js';

const state = {
  docs: [],
  active: 0,
  view: 'map',
  panel: 'Features',
  selection: null,
  highlight: null,
  enzymeSet: COMMON_SET,
  shownEnzymes: new Set(),
  digestWith: [],
  agarose: 1,
  ladder: '1 kb ladder',
  showTranslation: true,
  showComplement: true,
  showReads: true,
  trace: 0,
  pxPerBase: 14,
  cloning: {
    mode: 'restriction',
    vector: { docId: null, enzymes: [], fragment: null },
    insert: { docId: null, enzymes: [], fragment: null },
    circular: true,
    orientation: 'auto',
    gibson: [],
    minOverlap: 20,
  },
  restorable: null, // a saved session found on startup, until the user answers
};

const undo = new Map(); // doc id -> {past: [], future: []}
// Aligned reads are decoded from the file's own bytes rather than stored on the
// document, so they are never written to the session or saved back.
const alignments = new Map(); // doc id -> reads
const traces = new Map();     // doc id -> chromatograms

const readsFor = (doc) => (doc && state.showReads ? alignments.get(doc.id) ?? [] : []);

async function loadAlignment(doc) {
  if (alignments.has(doc.id)) return;
  alignments.set(doc.id, null);
  try {
    const reads = await readAlignments(doc);
    if (reads) { alignments.set(doc.id, reads); render(); }
  } catch (e) {
    toast(`Could not read the alignment in ${doc.name}: ${e.message}`, 'error');
  }
  try {
    const found = await readTraces(doc);
    if (found.length) { traces.set(doc.id, found); render(); }
  } catch (e) {
    toast(`Could not read the traces in ${doc.name}: ${e.message}`, 'error');
  }
}

// ---- remembering the session ------------------------------------------------

const isDirty = (doc) => !!(undo.get(doc.id)?.past.length || doc.restoredDirty);

let saveTimer = null;
function rememberSoon() {
  if (state.restorable) return; // nothing of the user's own is open yet
  clearTimeout(saveTimer);
  saveTimer = setTimeout(remember, 900);
}

function remember() {
  clearTimeout(saveTimer);
  if (!state.docs.length) { clearSession(); return; }
  saveSession({
    savedAt: new Date().toISOString(),
    active: state.active,
    view: state.view,
    panel: state.panel,
    selection: state.selection,
    enzymeSet: state.enzymeSet,
    shownEnzymes: [...state.shownEnzymes],
    digestWith: state.digestWith,
    agarose: state.agarose,
    ladder: state.ladder,
    showTranslation: state.showTranslation,
    showComplement: state.showComplement,
    // the undo history is deliberately left behind; only the documents return
    docs: state.docs.map((d) => ({ ...cloneDoc(d), restoredDirty: isDirty(d) })),
  });
}

function restoreSession(saved) {
  state.docs = saved.docs ?? [];
  for (const d of state.docs) undo.set(d.id, { past: [], future: [] });
  state.active = Math.min(saved.active ?? 0, Math.max(0, state.docs.length - 1));
  state.view = saved.view ?? 'map';
  state.panel = saved.panel ?? 'Features';
  state.selection = saved.selection ?? null;
  state.enzymeSet = saved.enzymeSet ?? COMMON_SET;
  state.shownEnzymes = new Set(saved.shownEnzymes ?? []);
  state.digestWith = saved.digestWith ?? [];
  state.agarose = saved.agarose ?? 1;
  state.ladder = saved.ladder ?? '1 kb ladder';
  state.showTranslation = saved.showTranslation ?? true;
  state.showComplement = saved.showComplement ?? true;
  state.restorable = null;
  render();
  for (const d of state.docs) loadAlignment(d);
  toast(`Reopened ${state.docs.length === 1 ? state.docs[0].name : `${state.docs.length} files`}`);
}

function startFresh() {
  state.restorable = null;
  clearSession();
  render();
}
const el = (id) => document.getElementById(id);
const current = () => state.docs[state.active] ?? null;

// ---- documents --------------------------------------------------------------

function openDoc(doc) {
  state.restorable = null; // opening something is an answer to the restore offer
  state.docs.push(doc);
  state.active = state.docs.length - 1;
  undo.set(doc.id, { past: [], future: [] });
  state.selection = null;
  state.highlight = null;
  state.digestWith = [];
  state.shownEnzymes = new Set();
  autoShowEnzymes(doc);
  render();
  loadAlignment(doc);
}

// Start by showing the unique cutters, which is what a map is usually for.
function autoShowEnzymes(doc) {
  const unique = [];
  for (const name of state.enzymeSet) {
    const n = findSites(doc.seq, doc.circular, ENZYME_BY_NAME[name]).length;
    if (n === 1) unique.push(name);
  }
  state.shownEnzymes = new Set(unique.slice(0, 14));
}

function closeDoc(i) {
  const doc = state.docs[i];
  if (doc && isDirty(doc) && !confirm(`Close ${doc.name}? Unsaved changes will be lost.`)) return;
  state.docs.splice(i, 1);
  undo.delete(doc?.id);
  state.active = Math.max(0, Math.min(state.active, state.docs.length - 1));
  state.selection = null;
  render();
}

function edit(label, mutate) {
  const doc = current();
  if (!doc) return;
  const stack = undo.get(doc.id);
  stack.past.push(cloneDoc(doc));
  if (stack.past.length > 60) stack.past.shift();
  stack.future.length = 0;
  mutate(doc);
  addHistory(doc, 'edit', label);
  render();
}

function undoOnce(redo = false) {
  const doc = current();
  if (!doc) return;
  const stack = undo.get(doc.id);
  const from = redo ? stack.future : stack.past;
  if (!from.length) return;
  const to = redo ? stack.past : stack.future;
  to.push(cloneDoc(doc));
  const prev = from.pop();
  state.docs[state.active] = prev;
  undo.set(prev.id, stack);
  state.selection = null;
  render();
}

// ---- files ------------------------------------------------------------------

async function openFiles(files) {
  for (const file of files) {
    try {
      const lower = file.name.toLowerCase();
      if (lower.endsWith('.dna')) {
        const doc = readDna(new Uint8Array(await file.arrayBuffer()), file.name);
        openDoc(doc);
        for (const w of doc.source.warnings ?? []) toast(w);
      } else {
        const text = await file.text();
        const docs = /^LOCUS/m.test(text) ? readGenBank(text, file.name)
          : text.trimStart().startsWith('>') ? readFasta(text, file.name)
            : readFasta(`>${file.name}\n${text}`, file.name);
        for (const d of docs) openDoc(d);
      }
    } catch (e) {
      toast(`${file.name}: ${e.message}`, 'error');
    }
  }
}

function save(kind) {
  const doc = current();
  if (!doc) return;
  let data;
  let name = doc.name.replace(/[^\w.-]+/g, '_');
  if (kind === 'dna') {
    data = new Blob([writeDna(doc, { bindingSitesFor: (p) => findBindingSites(doc.seq, doc.circular, p.seq) })], { type: 'application/octet-stream' });
    name += '.dna';
  } else if (kind === 'gb') {
    data = new Blob([writeGenBank(doc)], { type: 'chemical/seq-na-genbank' });
    name += '.gb';
  } else {
    data = new Blob([writeFasta(doc)], { type: 'text/plain' });
    name += '.fasta';
  }
  const url = URL.createObjectURL(data);
  const a = h('a', { href: url, download: name });
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 4000);
  toast(`Saved ${name}`);
}

// ---- api handed to the panels ------------------------------------------------

const api = {
  edit,
  openDoc,
  setSelection(sel) { state.selection = sel; render(); },
  focusFeature(f) {
    state.highlight = state.highlight === f.id ? null : f.id;
    if (state.highlight) {
      state.selection = { start: featureSpan(f).start, end: featureSpan(f).end };
      state.panel = 'Features'; // so the details appear in the sidebar
    }
    render();
  },
  toggleEnzyme(name) {
    state.shownEnzymes.has(name) ? state.shownEnzymes.delete(name) : state.shownEnzymes.add(name);
    render();
  },
  showEnzymes(names) { state.shownEnzymes = new Set(names); render(); },
  setEnzymeSet(names) { state.enzymeSet = names; render(); },
  toggleDigest(name) {
    state.digestWith = state.digestWith.includes(name) ? state.digestWith.filter((x) => x !== name) : [...state.digestWith, name];
    render();
  },
  setAgarose(v) { state.agarose = v; render(); },
  setLadder(v) { state.ladder = v; render(); },
  primerFromSelection(dir) { addPrimerFromSelection(dir); },
  setCloning(patch) { state.cloning = { ...state.cloning, ...patch }; render(); },
  selectedSequence() { return selectedSeq(); },
  copy(text, what) { copy(text, what); },
};

// ---- rendering ---------------------------------------------------------------

function shownSites(doc) {
  const out = [];
  for (const name of state.shownEnzymes) {
    const e = ENZYME_BY_NAME[name];
    if (e) out.push(...findSites(doc.seq, doc.circular, e));
  }
  return out.sort((a, b) => a.top - b.top);
}

function renderTabs() {
  const box = clear(el('tabs'));
  state.docs.forEach((doc, i) => {
    const dirty = isDirty(doc);
    box.append(h('div', { class: 'tab', role: 'tab', 'aria-selected': String(i === state.active), onclick: () => { state.active = i; state.selection = null; state.highlight = null; render(); } },
      h('span', { class: 'name' }, doc.name),
      dirty ? h('span', { class: 'dirty', title: 'edited' }, '•') : null,
      h('button', { class: 'close', title: 'Close', onclick: (e) => { e.stopPropagation(); closeDoc(i); } }, '×')));
  });
}

function renderViewBar(doc) {
  const opts = clear(el('viewOpts'));
  const ops = clear(el('editOps'));
  if (!doc) return;
  const toggles = [['showComplement', 'Both strands'], ['showTranslation', 'Translation']];
  if (alignments.get(doc.id)?.length) toggles.push(['showReads', `Reads (${alignments.get(doc.id).length})`]);
  for (const [key, label] of toggles) {
    opts.append(h('label', { class: 'sub', style: { display: 'flex', gap: '4px', alignItems: 'center' } },
      h('input', { type: 'checkbox', checked: state[key], onchange: (e) => { state[key] = e.target.checked; render(); } }), label));
  }
  const sel = state.selection;
  ops.append(h('button', {
    class: 'btn', disabled: !sel,
    onclick: () => edit('Delete selection', (d) => { replaceRange(d, sel.start, sel.end, ''); }),
  }, 'Delete'));
  ops.append(h('button', {
    class: 'btn', disabled: !sel,
    onclick: () => {
      const seq = cleanSequence(prompt?.('Replace the selection with:') ?? '');
      edit('Replace selection', (d) => { replaceRange(d, sel.start, sel.end, seq); });
    },
  }, 'Replace…'));
  ops.append(h('button', {
    class: 'btn', disabled: !sel || !doc.circular,
    title: 'Make the selection start base 1',
    onclick: () => edit('Set origin', (d) => { rotateOrigin(d, sel.start); }),
  }, 'Set origin'));
  ops.append(h('button', { class: 'btn', onclick: () => edit('Reverse complement', (d) => reverseComplementDoc(d)) }, 'Rev. comp.'));
}

function renderStatus(doc) {
  const bar = clear(el('status'));
  if (!doc) return;
  bar.append(h('span', {}, h('b', {}, fmt(doc.seq.length)), ' bp · ', doc.circular ? 'circular' : 'linear'));
  bar.append(h('span', {}, 'GC ', h('b', {}, `${gcContent(doc.seq).toFixed(1)}%`)));
  const sel = state.selection;
  if (sel) {
    const len = rangeLength(sel, doc.seq.length);
    const s = subseq(doc.seq, sel.start, sel.end, doc.circular);
    bar.append(h('span', {}, 'Selection ', h('b', {}, `${fmt(sel.start + 1)}–${fmt(sel.end)}`), ` · ${fmt(len)} bp · GC ${gcContent(s).toFixed(1)}%`));
    if (len >= 8 && len <= 60) bar.append(h('span', {}, 'Tm ', h('b', {}, `${meltingTemp(s).toFixed(1)}°C`)));
    bar.append(h('button', { class: 'btn', style: { padding: '2px 8px' }, onclick: () => copy(s) }, 'Copy'));
    bar.append(h('button', { class: 'btn', style: { padding: '2px 8px' }, onclick: () => copy(revcomp(s)) }, 'Copy rev. comp.'));
  } else {
    bar.append(h('span', {}, 'Drag to select · right-click for everything you can do with it'));
  }
}

// Copy through a scratch textarea first: it is synchronous, so it still counts
// as part of the key press or click that asked for it. The async clipboard API
// is the fallback, since some browsers refuse it outside a gesture.
function copyViaTextarea(text) {
  try {
    const ta = document.createElement('textarea');
    ta.value = text;
    ta.setAttribute('readonly', '');
    ta.style.cssText = 'position:fixed;top:0;left:-9999px;opacity:0';
    document.body.append(ta);
    ta.select();
    const ok = document.execCommand('copy');
    ta.remove();
    return ok;
  } catch {
    return false;
  }
}

function copy(text, what = 'bases') {
  if (!text) return;
  const done = `Copied ${fmt(text.length)} ${what}`;
  if (copyViaTextarea(text)) { toast(done); return; }
  navigator.clipboard?.writeText(text).then(
    () => toast(done),
    () => toast('Could not reach the clipboard', 'error'),
  ) ?? toast('Could not reach the clipboard', 'error');
}

function selectedSeq() {
  const doc = current();
  if (!doc || !state.selection) return '';
  return subseq(doc.seq, state.selection.start, state.selection.end, doc.circular);
}

function addPrimerFromSelection(dir) {
  const picked = selectedSeq();
  const sel = state.selection;
  if (!picked) return;
  state.panel = 'Primers';
  edit(`Add ${dir === 1 ? 'forward' : 'reverse'} primer`, (d) => {
    d.primers.push({
      id: newId(),
      name: `${dir === 1 ? 'fwd' : 'rev'} ${sel.start + 1}-${sel.end}`,
      seq: dir === 1 ? picked : revcomp(picked),
      description: `From ${sel.start + 1}–${sel.end} of ${d.name}`,
      extra: {},
    });
  });
}

function addFeatureFromSelection() {
  const sel = state.selection;
  if (!sel) return;
  const id = newId();
  state.panel = 'Features';
  state.highlight = id;
  edit('Add feature', (d) => {
    d.features.push({
      id, name: 'new feature', type: 'misc_feature', strand: 1, color: colorFor('misc_feature'),
      segments: [{ start: sel.start, end: sel.end }], qualifiers: [], extra: {},
    });
  });
}

function featureFromEvent(ev) {
  const id = ev.target?.closest?.('[data-feature]')?.getAttribute('data-feature');
  return id ? current()?.features.find((f) => f.id === id) : null;
}

function openContextMenu(ev) {
  const doc = current();
  if (!doc) return;
  ev.preventDefault();
  const feature = featureFromEvent(ev);
  const sel = state.selection;
  const picked = selectedSeq();
  const items = [];

  if (feature) {
    items.push({ label: `Select “${feature.name}”`, onClick: () => api.focusFeature(feature) });
    items.push({ label: 'Copy feature sequence', onClick: () => copy(featureSequence(doc, feature)) });
    if (feature.type === 'CDS') {
      items.push({ label: 'Copy protein', onClick: () => copy(translate(codingSequence(doc, feature)), 'residues') });
    }
    items.push({ label: 'Delete feature', danger: true, onClick: () => edit(`Delete ${feature.name}`, (d) => { d.features = d.features.filter((x) => x.id !== feature.id); }) });
    items.push('-');
  }

  items.push({ label: 'Copy', hint: '⌘C', disabled: !picked, onClick: () => copy(picked) });
  items.push({ label: 'Copy reverse complement', disabled: !picked, onClick: () => copy(revcomp(picked)) });
  items.push({ label: 'Copy translation', disabled: !picked, onClick: () => copy(translate(picked), 'residues') });
  items.push({
    label: 'Copy as FASTA', disabled: !picked,
    onClick: () => copy(`>${doc.name} ${sel.start + 1}-${sel.end}\n${picked.match(/.{1,70}/g).join('\n')}\n`),
  });
  items.push('-');
  items.push({ label: 'Make a forward primer', disabled: !picked, onClick: () => addPrimerFromSelection(1) });
  items.push({ label: 'Make a reverse primer', disabled: !picked, onClick: () => addPrimerFromSelection(-1) });
  items.push({ label: 'Design a primer pair…', disabled: !picked, onClick: () => { state.panel = 'Primers'; render(); } });
  items.push('-');
  items.push({ label: 'Add a feature here', disabled: !sel, onClick: addFeatureFromSelection });
  items.push({ label: 'Delete selected bases', danger: true, disabled: !sel, onClick: () => edit('Delete selection', (d) => replaceRange(d, sel.start, sel.end, '')) });
  items.push({
    label: 'Replace selected bases…', disabled: !sel,
    onClick: () => {
      const text = cleanSequence(window.prompt('Replace the selection with:') ?? '');
      edit('Replace selection', (d) => replaceRange(d, sel.start, sel.end, text));
    },
  });
  items.push('-');
  if (doc.circular && sel) items.push({ label: `Start the sequence at ${fmt(sel.start + 1)}`, onClick: () => edit('Set origin', (d) => rotateOrigin(d, sel.start)) });
  items.push({ label: 'Select everything', hint: '⌘A', onClick: () => api.setSelection({ start: 0, end: doc.seq.length }) });
  items.push({ label: 'Clear the selection', disabled: !sel, onClick: () => api.setSelection(null) });
  showMenu(ev.clientX, ev.clientY, items);
}

function renderView(doc) {
  // keep the reader where they were: rebuilding the view resets its scroll
  const wasScrolled = state.view === 'seq' ? el('view').scrollTop : 0;
  const view = clear(el('view'));
  if (!doc) { view.append(welcome()); return; }
  const sites = shownSites(doc);
  const opts = {
    selection: state.selection, sites, highlight: state.highlight,
    onPick: (f) => f && api.focusFeature(f),
    onPickSite: (st) => api.setSelection({ start: st.site.start, end: st.site.end }),
  };
  if (state.view === 'trace') {
    const found = traces.get(doc.id) ?? [];
    const tv = traceView(found, {
      selected: Math.min(state.trace, found.length - 1),
      pxPerBase: state.pxPerBase,
      onSelect: (i) => { state.trace = i; render(); },
      onZoom: (px) => { state.pxPerBase = px; render(); },
      onPickBase: (k, trace) => {
        const base = trace.bases[k];
        toast(`Base ${k + 1} of ${trace.name || 'the trace'}: ${base}${trace.confidence ? ` (quality ${trace.confidence[k] ?? '?'})` : ''}`);
      },
    });
    view.append(tv);
    addEventListener('resize', () => tv.refresh?.(), { once: true });
    return;
  }
  if (state.view === 'map') {
    const wrap = h('div', { class: 'mapwrap' });
    wrap.append(doc.circular ? circularMap(doc, opts) : linearMap(doc, opts));
    wrap.addEventListener('click', (e) => { if (e.target === wrap) { state.highlight = null; state.selection = null; render(); } });
    view.append(wrap);
  } else {
    const sv = sequenceView(doc, {
      ...opts,
      // While the pointer is down the view repaints itself; rebuilding it here
      // would drop the drag. Only once the drag ends do the panels catch up.
      onSelect: (sel, info) => {
        state.selection = sel;
        renderStatus(doc);
        if (!info?.live) render();
      },
      showTranslation: state.showTranslation, showComplement: state.showComplement,
      reads: readsFor(doc),
    });
    view.append(sv);
    view.onscroll = () => sv.onScroll?.();
    if (wasScrolled) {
      requestAnimationFrame(() => { view.scrollTop = wasScrolled; sv.onScroll?.(); });
    }
  }
}

function renderPanels(doc) {
  const tabs = clear(el('panelTabs'));
  for (const name of PANELS) {
    tabs.append(h('button', {
      role: 'tab', 'aria-selected': String(state.panel === name),
      onclick: () => { state.panel = name; render(); },
    }, name));
  }
  const panel = clear(el('panel'));
  if (!doc) { panel.append(h('p', { class: 'empty' }, 'Open a file to see its details.')); return; }
  panel.append(renderPanel(state.panel, { doc, state, api, docs: state.docs }));
}

function welcome() {
  if (state.restorable) return restoreCard(state.restorable);
  const drop = h('div', { class: 'drop' },
    h('p', { style: { margin: '0 0 12px' } }, 'Drop a sequence file here'),
    h('button', { class: 'btn primary', onclick: () => el('fileInput').click() }, 'Choose a file…'));
  const box = h('div', { class: 'welcome' },
    h('h1', {}, 'A plasmid map in the browser'),
    h('p', {}, 'Open a SnapGene .dna, GenBank or FASTA file to see its map, sequence, features and cut sites. Files stay on this computer.'),
    drop,
    h('div', { class: 'group' },
      h('button', { class: 'btn', onclick: () => openDoc(exampleDocument()) }, 'Open the example fragment'),
      h('button', { class: 'btn', onclick: () => openDoc(newDocument({ name: 'Untitled', seq: '', circular: true })) }, 'Start an empty sequence')),
    h('ul', { style: { marginTop: '22px' } },
      h('li', {}, 'Click a feature or an enzyme label to select it'),
      h('li', {}, 'Drag across the sequence to select a region, then ⌘C to copy it'),
      h('li', {}, 'Right-click anywhere for copying, primers and editing'),
      h('li', {}, 'Saving a .dna file keeps the parts BublaGene does not edit, such as the cloning history')),
    h('p', { class: 'sub', style: { marginTop: '18px' } },
      'Free software under the AGPL-3.0. ',
      h('a', { href: 'https://github.com/KowalskiBio/bublagene', target: '_blank', rel: 'noopener', style: { color: 'var(--accent)' } }, 'Source code')));
  return box;
}

function restoreCard(saved) {
  const docs = saved.docs ?? [];
  const list = h('ul', { style: { listStyle: 'none', padding: '0', margin: '0 0 20px' } });
  for (const d of docs) {
    list.append(h('li', {
      style: { display: 'flex', gap: '10px', alignItems: 'baseline', padding: '7px 0', borderBottom: '1px solid var(--line)' },
    },
      h('span', { style: { fontWeight: '500' } }, d.name),
      h('span', { class: 'sub' }, `${fmt((d.seq ?? '').length)} bp · ${d.circular ? 'circular' : 'linear'}`),
      d.restoredDirty ? h('span', { class: 'sub', style: { color: 'var(--accent)' } }, '· unsaved changes') : null));
  }
  return h('div', { class: 'welcome' },
    h('h1', {}, 'Pick up where you left off?'),
    h('p', {}, `${docs.length === 1 ? 'One file was' : `${docs.length} files were`} open ${timeAgo(saved.savedAt)}.`),
    list,
    h('div', { class: 'group' },
      h('button', { class: 'btn primary', onclick: () => restoreSession(saved) }, docs.length === 1 ? 'Reopen it' : 'Reopen them'),
      h('button', { class: 'btn', onclick: startFresh }, 'Start fresh')),
    h('p', { class: 'sub', style: { marginTop: '18px' } },
      'Kept in this browser only, never uploaded. “Start fresh” forgets it.'));
}

function render() {
  const doc = current();
  renderTabs();
  const hasTraces = !!traces.get(doc?.id)?.length;
  for (const b of el('viewSeg').querySelectorAll('button')) {
    b.setAttribute('aria-pressed', String(b.dataset.view === state.view));
    if (b.dataset.view === 'trace') b.hidden = !hasTraces;
  }
  if (state.view === 'trace' && !hasTraces) state.view = 'map';
  renderViewBar(doc);
  renderView(doc);
  renderStatus(doc);
  renderPanels(doc);
  el('undoBtn').disabled = !(undo.get(doc?.id)?.past.length);
  el('redoBtn').disabled = !(undo.get(doc?.id)?.future.length);
  el('saveAs').disabled = !doc;
  document.title = doc ? `${doc.name} · BublaGene` : 'BublaGene';
  rememberSoon();
}

// ---- events -------------------------------------------------------------------

el('openBtn').addEventListener('click', () => el('fileInput').click());
el('fileInput').addEventListener('change', (e) => { openFiles([...e.target.files]); e.target.value = ''; });
el('newBtn').addEventListener('click', () => openDoc(newDocument({ name: 'Untitled', seq: '', circular: true })));
el('saveAs').addEventListener('change', (e) => { if (e.target.value) save(e.target.value); e.target.value = ''; });
el('undoBtn').addEventListener('click', () => undoOnce(false));
el('redoBtn').addEventListener('click', () => undoOnce(true));
el('themeBtn').addEventListener('click', () => {
  const now = document.documentElement.dataset.theme;
  const next = now === 'dark' ? 'light' : now === 'light' ? '' : 'dark';
  if (next) document.documentElement.dataset.theme = next;
  else delete document.documentElement.dataset.theme;
  try { localStorage.setItem('bublagene.theme', next); } catch { /* ignore */ }
});
el('viewSeg').addEventListener('click', (e) => {
  const v = e.target.dataset?.view;
  if (v) { state.view = v; render(); }
});

document.addEventListener('dragover', (e) => { e.preventDefault(); document.querySelector('.drop')?.classList.add('over'); });
document.addEventListener('dragleave', () => document.querySelector('.drop')?.classList.remove('over'));
document.addEventListener('drop', (e) => {
  e.preventDefault();
  document.querySelector('.drop')?.classList.remove('over');
  if (e.dataTransfer?.files?.length) openFiles([...e.dataTransfer.files]);
});

el('view').addEventListener('contextmenu', openContextMenu);

document.addEventListener('keydown', (e) => {
  const typing = /^(INPUT|TEXTAREA|SELECT)$/.test(e.target.tagName) || e.target.isContentEditable;
  const mod = e.metaKey || e.ctrlKey;
  const key = e.key.toLowerCase();
  if (mod && key === 'z') { e.preventDefault(); undoOnce(e.shiftKey); }
  else if (mod && key === 'o') { e.preventDefault(); el('fileInput').click(); }
  else if (mod && key === 'c' && !typing && state.selection) { e.preventDefault(); copy(selectedSeq()); }
  else if (mod && key === 'a' && !typing && current()) { e.preventDefault(); api.setSelection({ start: 0, end: current().seq.length }); }
  else if (!typing && e.key === 'Escape') { closeMenu(); state.selection = null; state.highlight = null; render(); }
});

// the browser's own Copy (menu bar, or a trackpad gesture) — our selection is
// not a native text selection, so hand it the bases ourselves
document.addEventListener('copy', (e) => {
  if (/^(INPUT|TEXTAREA)$/.test(e.target.tagName) || e.target.isContentEditable) return;
  const picked = selectedSeq();
  if (!picked) return;
  e.preventDefault();
  e.clipboardData.setData('text/plain', picked);
  toast(`Copied ${fmt(picked.length)} bases`);
});

window.addEventListener('resize', () => { if (state.view === 'seq') renderView(current()); });

// write the session out before the tab goes away, rather than waiting out the debounce
addEventListener('pagehide', remember);
document.addEventListener('visibilitychange', () => { if (document.visibilityState === 'hidden') remember(); });

try {
  const t = localStorage.getItem('bublagene.theme');
  if (t) document.documentElement.dataset.theme = t;
} catch { /* ignore */ }

render();

// offer the previous session once storage answers
loadSession().then((saved) => {
  if (saved?.docs?.length && !state.docs.length) {
    state.restorable = saved;
    render();
  }
});
