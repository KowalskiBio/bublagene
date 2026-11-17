// @ts-nocheck — legacy DOM wiring; the engines it calls are type-checked.
// BublaGene: application state and wiring.

import { h, clear, fmt, toast } from './dom.js';
import { circularMap, linearMap } from './map.js';
import { sequenceView } from './seqview.js';
import { showMenu, closeMenu } from './menu.js';
import { saveSession, loadSession, clearSession, timeAgo } from './session.js';
import { readAlignments } from '../formats/bam.js';
import { readTraces } from '../formats/ztr.js';
import { traceBaseToReference } from '../traces.js';
import { traceView } from './traceview.js';
import { History } from '../undo.js';
import { PANELS, renderPanel } from './panels.js';
import { exampleDocument } from './sample.js';
import { readDna, writeDna } from '../formats/dna.js';
import { readGenBank, writeGenBank, readFasta, writeFasta } from '../formats/genbank.js';
import { newDocument, cloneDoc, replaceRange, rotateOrigin, reverseComplementDoc, addHistory, featureSpan, featureSequence, codingSequence, newId, colorFor } from '../model.js';
import { findSites, ENZYME_BY_NAME, COMMON_SET } from '../enzymes.js';
import { findBindingSites } from '../primers.js';
import { gcContent, meltingTemp, revcomp, subseq, rangeLength, cleanSequence, translate } from '../seq.js';
import { alignPair } from '../align.js';
import { callVariants, variantsToFeatures } from '../variants.js';
import { findGuides, guidesToFeatures, nhejEdit, hdrEdit, editPrimers } from '../crispr.js';
import { fold, structureToSvg } from '../fold.js';
import { extractFragment, pasteFragment, reorderFeature } from '../fragment.js';
import { offload } from '../compute.js';
import { saveDocument, openDocument, supportsFileSystemAccess, writeAutosave, readAutosave } from '../storage.js';
import { projectToJSON, projectFromJSON } from '../project.js';

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
  showSixFrame: false,
  showComplement: true,
  featureFilter: '',
  showReads: true,
  trace: 0,
  traceBase: null,
  selectedFeatures: new Set(),
  alignTarget: null,
  alignLocal: false,
  alignBand: 0,
  alignBusy: false,
  alignResult: null,
  crispr: null,
  fold: null,
  clipboard: null,
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
const traceMaps = new Map();   // doc id -> { key, map }
let activeTraceView = null;

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

const isDirty = (doc) => !!(undo.get(doc.id)?.canUndo?.() || doc.restoredDirty);

let saveTimer = null;
function rememberSoon() {
  if (state.restorable) return; // nothing of the user's own is open yet
  clearTimeout(saveTimer);
  saveTimer = setTimeout(remember, 900);
}

function remember() {
  clearTimeout(saveTimer);
  if (!state.docs.length) { clearSession(); return; }
  const payload = {
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
    showSixFrame: state.showSixFrame,
    showComplement: state.showComplement,
    featureFilter: state.featureFilter,
    // the undo history is deliberately left behind; only the documents return
    docs: state.docs.map((d) => ({ ...cloneDoc(d), restoredDirty: isDirty(d) })),
  };
  saveSession(payload);
  // A small JSON mirror for recovery when IndexedDB is unavailable or wedged.
  try {
    const backup = {
      ...payload,
      docs: payload.docs.map((d) => ({ ...d, source: { format: 'backup' } })),
    };
    writeAutosave('bublagene.backup', JSON.stringify(backup));
  } catch { /* storage full or unavailable */ }
}

function restoreSession(saved) {
  state.docs = saved.docs ?? [];
  for (const d of state.docs) undo.set(d.id, new History());
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
  state.showSixFrame = saved.showSixFrame ?? false;
  state.showComplement = saved.showComplement ?? true;
  state.featureFilter = saved.featureFilter ?? '';
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
  undo.set(doc.id, new History());
  state.selection = null;
  state.highlight = null;
  state.selectedFeatures = new Set();
  state.traceBase = null;
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
  const wasActive = i === state.active;
  state.docs.splice(i, 1);
  undo.delete(doc?.id);
  // closing a tab to the left shifts every later tab down; keep the same document open
  if (i < state.active) state.active -= 1;
  state.active = Math.max(0, Math.min(state.active, state.docs.length - 1));
  if (wasActive) {
    state.selection = null;
    state.highlight = null;
    state.selectedFeatures = new Set();
    state.traceBase = null;
  }
  render();
}

function edit(label, mutate) {
  const doc = current();
  if (!doc) return;
  const stack = undo.get(doc.id);
  stack.record(cloneDoc(doc), label, { coalesceKey: label });
  mutate(doc);
  state.alignResult = null;
  state.fold = null;
  state.crispr = null;
  // a selection that no longer fits the edited sequence would point past its end
  if (state.selection && (state.selection.end > doc.seq.length || state.selection.start > doc.seq.length)) state.selection = null;
  addHistory(doc, 'edit', label);
  render();
}

// Replace the selection with typed text. Cancelling the prompt must leave the
// sequence alone; an empty answer is treated the same way (Delete is for that).
function replaceSelection() {
  const sel = state.selection;
  if (!sel) return;
  const answer = window.prompt('Replace the selection with:');
  if (answer === null) return;
  const seq = cleanSequence(answer);
  if (!seq) return;
  edit('Replace selection', (d) => { replaceRange(d, sel.start, sel.end, seq); });
  state.selection = { start: sel.start, end: sel.start + seq.length };
  render();
}

function deleteSelection() {
  const sel = state.selection;
  if (!sel) return;
  edit('Delete selection', (d) => { replaceRange(d, sel.start, sel.end, ''); });
  state.selection = null;
  render();
}

// A highlighted feature is only "current" while the selection is that feature;
// once the reader selects something else, or the feature is gone, drop it.
function syncHighlight() {
  if (!state.highlight) return;
  const f = current()?.features.find((x) => x.id === state.highlight);
  const sel = state.selection;
  const span = f && featureSpan(f);
  if (!f || !sel || sel.start !== span.start || sel.end !== span.end) state.highlight = null;
}

function undoOnce(redo = false) {
  const doc = current();
  if (!doc) return;
  const stack = undo.get(doc.id);
  const prev = redo ? stack.redo(cloneDoc(doc)) : stack.undo(cloneDoc(doc));
  if (!prev) return;
  state.docs[state.active] = prev;
  undo.set(prev.id, stack);
  state.selection = null;
  state.selectedFeatures = new Set();
  render();
}

// ---- files ------------------------------------------------------------------

async function readFileData(file) {
  if (file.data) return file.data instanceof Uint8Array ? file.data : new Uint8Array(file.data);
  return new Uint8Array(await file.arrayBuffer());
}

async function readFileText(file) {
  if (file.text) return file.text();
  return new TextDecoder().decode(await readFileData(file));
}

async function openFiles(files) {
  for (const file of files) {
    try {
      const lower = file.name.toLowerCase();
      if (lower.endsWith('.dna')) {
        const doc = readDna(await readFileData(file), file.name);
        openDoc(doc);
        for (const w of doc.source.warnings ?? []) toast(w);
      } else if (lower.endsWith('.json')) {
        const project = projectFromJSON(await readFileText(file));
        for (const d of project.docs) openDoc(d);
        if (project.warnings.length) toast(project.warnings[0]);
      } else {
        const text = await readFileText(file);
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

async function save(kind) {
  const doc = current();
  if (!doc) return;
  let data;
  let name = doc.name.replace(/[^\w.-]+/g, '_');
  let types;
  if (kind === 'dna') {
    data = writeDna(doc, { bindingSitesFor: (p) => findBindingSites(doc.seq, doc.circular, p.seq) });
    name += '.dna';
    types = [{ description: 'BublaGene sequence', accept: { 'application/octet-stream': ['.dna'] } }];
  } else if (kind === 'gb') {
    data = writeGenBank(doc);
    name += '.gb';
    types = [{ description: 'GenBank', accept: { 'chemical/seq-na-genbank': ['.gb', '.gbk'] } }];
  } else {
    data = writeFasta(doc);
    name += '.fasta';
    types = [{ description: 'FASTA', accept: { 'text/plain': ['.fasta', '.fa'] } }];
  }
  try {
    const result = await saveDocument({ data, name, types });
    toast(result.method === 'download' ? `Saved ${result.name}` : `${result.name} saved`);
  } catch (e) {
    toast(`Could not save: ${e.message}`, 'error');
  }
}

async function saveProject() {
  if (!state.docs.length) return;
  const json = projectToJSON(state.docs, {
    active: state.active,
    settings: {
      view: state.view, panel: state.panel, enzymeSet: state.enzymeSet,
      shownEnzymes: [...state.shownEnzymes], digestWith: state.digestWith,
      agarose: state.agarose, ladder: state.ladder,
      showTranslation: state.showTranslation, showSixFrame: state.showSixFrame, showComplement: state.showComplement,
      featureFilter: state.featureFilter,
    },
    alignments: new Map([...alignments.entries()].filter(([, reads]) => reads).map(([id, reads]) => [id, reads.map((r) => ({ ...r, insertions: r.insertions instanceof Map ? [...r.insertions.entries()] : (r.insertions ?? []) }))])),
  });
  try {
    const result = await saveDocument({ data: json, name: 'project.json', types: [{ description: 'BublaGene project', accept: { 'application/json': ['.json'] } }] });
    toast(result.method === 'download' ? 'Project downloaded' : 'Project saved');
  } catch (e) { toast(`Could not save the project: ${e.message}`, 'error'); }
}

async function openProject() {
  try {
    const file = await openDocument({ types: [{ description: 'BublaGene project', accept: { 'application/json': ['.json'] } }] });
    if (!file) return;
    const project = projectFromJSON(new TextDecoder().decode(file.data));
    for (const entry of project.alignments ?? []) {
      alignments.set(entry.docId, (entry.reads ?? []).map((r) => ({ ...r, insertions: r.insertions instanceof Map ? r.insertions : new Map(r.insertions ?? []) })));
    }
    for (const d of project.docs) openDoc(d);
    const s = project.settings ?? {};
    state.view = s.view ?? state.view;
    state.panel = s.panel ?? state.panel;
    state.enzymeSet = s.enzymeSet ?? state.enzymeSet;
    state.shownEnzymes = new Set(s.shownEnzymes ?? state.shownEnzymes);
    state.digestWith = s.digestWith ?? state.digestWith;
    state.showTranslation = s.showTranslation ?? state.showTranslation;
    state.showSixFrame = s.showSixFrame ?? state.showSixFrame;
    state.showComplement = s.showComplement ?? state.showComplement;
    for (const w of project.warnings) toast(w);
    toast(`Opened ${project.docs.length} document(s)`);
    render();
  } catch (e) { toast(`Could not open the project: ${e.message}`, 'error'); }
}

// ---- api handed to the panels ------------------------------------------------

const api = {
  edit,
  openDoc,
  setSelection(sel) { state.selection = sel; syncHighlight(); render(); },
  setFeatureFilter(text) { state.featureFilter = text; render(); },
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
  // quiet: the caller redraws the gel itself, so a slider being dragged is not rebuilt under the pointer
  setAgarose(v, quiet = false) { state.agarose = v; if (quiet) rememberSoon(); else render(); },
  clearCrispr() { state.crispr = null; render(); },
  setLadder(v) { state.ladder = v; render(); },
  primerFromSelection(dir) { addPrimerFromSelection(dir); },
  setCloning(patch) { state.cloning = { ...state.cloning, ...patch }; render(); },
  selectedSequence() { return selectedSeq(); },
  copy(text, what) { copy(text, what); },
  featureSequence(f) { const doc = current(); return doc ? featureSequence(doc, f) : ''; },
  toggleFeatureSelected(id) {
    state.selectedFeatures.has(id) ? state.selectedFeatures.delete(id) : state.selectedFeatures.add(id);
    render();
  },
  setFeatureSelection(ids) { state.selectedFeatures = new Set(ids); render(); },
  deleteSelectedFeatures() {
    const doc = current();
    if (!doc || !state.selectedFeatures.size) return;
    const ids = new Set(state.selectedFeatures);
    edit('Delete selected features', (d) => { d.features = d.features.filter((f) => !ids.has(f.id)); });
    state.selectedFeatures = new Set();
    render();
  },
  copySelectedFeatures() {
    const doc = current();
    if (!doc) return;
    const picked = doc.features.filter((f) => state.selectedFeatures.has(f.id));
    copy(picked.map((f) => `>${f.name}\n${featureSequence(doc, f)}`).join('\n'), 'feature sequences');
  },
  reorderFeature(id, beforeId) { edit('Reorder features', (d) => reorderFeature(d, id, beforeId)); },
  copyFragment() {
    const doc = current();
    if (!doc || !state.selection) return;
    state.clipboard = extractFragment(doc, state.selection);
    if (state.clipboard) toast(`Copied ${fmt(state.clipboard.seq.length)} bp with ${state.clipboard.features.length} feature(s)`);
  },
  pasteFragment() {
    const doc = current();
    if (!doc || !state.clipboard) return;
    const at = state.selection?.start ?? 0;
    edit('Paste fragment', (d) => pasteFragment(d, at, state.clipboard));
    toast(`Pasted ${fmt(state.clipboard.seq.length)} bp`);
  },
  setAlignTarget(id) { state.alignTarget = id; state.alignResult = null; render(); },
  setAlignLocal(v) { state.alignLocal = v; render(); },
  setAlignBand(v) { state.alignBand = v; render(); },
  async runAlignment() {
    const doc = current();
    const others = state.docs.filter((d) => d.id !== doc?.id);
    const target = others.find((d) => d.id === state.alignTarget) ?? others[0];
    if (!doc || !target) return;
    state.alignBusy = true;
    render();
    try {
      const options = { local: state.alignLocal, band: state.alignBand > 0 ? state.alignBand : null };
      const { result } = await offload(
        ({ a, b, options: opts }) => alignPair(a, b, opts),
        { kind: 'align', a: doc.seq, b: target.seq, options },
        { workerUrl: new URL('../workers/compute-worker.js', import.meta.url) },
      );
      if (result?.__error) throw new Error(result.__error);
      state.alignResult = { docId: doc.id, targetId: target.id, result };
    } catch (e) {
      state.alignResult = { docId: doc.id, targetId: target.id, error: e.message };
    }
    state.alignBusy = false;
    render();
  },
  alignmentVariants() {
    const doc = current();
    const result = state.alignResult;
    const target = state.docs.find((d) => d.id === result?.targetId);
    if (!doc || !target) return;
    const calls = callVariants(target.seq, doc.seq);
    edit('Add variant features', (d) => { d.features.push(...variantsToFeatures(d, calls.variants)); });
    toast(`${calls.variants.length} variant(s) added as features`);
  },
  runCrispr(opts) {
    const doc = current();
    if (!doc) return;
    try {
      const guides = findGuides(doc.seq, {
        enzyme: opts.enzyme, guideLength: opts.guideLength, minGC: opts.minGC, maxGC: opts.maxGC,
        circular: doc.circular, maxOffTargets: opts.offTargets ? 3 : null,
      });
      state.crispr = { docId: doc.id, ...opts, guides, error: null };
    } catch (e) {
      state.crispr = { docId: doc.id, ...opts, guides: [], error: e.message };
    }
    render();
  },
  addCrisprFeatures() {
    const doc = current();
    const guides = state.crispr?.docId === doc?.id ? state.crispr.guides : [];
    if (!guides.length) return;
    edit('Add CRISPR guide features', (d) => { d.features.push(...guidesToFeatures(d, guides)); });
  },
  simulateNhej(guide) {
    const doc = current();
    if (!doc) return;
    const result = nhejEdit(doc, guide, { seed: 1 });
    openDoc(result.doc);
    toast(`NHEJ: deleted ${result.change.deleted.length} bp`);
  },
  simulateHdr(guide, donor) {
    const doc = current();
    if (!doc) return;
    const result = hdrEdit(doc, guide, { donor });
    openDoc(result.doc);
    toast(`HDR: inserted ${donor.length} bp`);
  },
  crisprPrimers(guide) {
    const doc = current();
    if (!doc) return;
    try {
      const primers = editPrimers(doc, guide.cut, { arm: 40 });
      edit('Add edit primers', (d) => {
        d.primers.push({ id: newId(), name: 'edit forward', seq: primers.forward.sequence, description: `Flanks CRISPR cut at ${guide.cut + 1}`, extra: {} });
        d.primers.push({ id: newId(), name: 'edit reverse', seq: primers.reverse.sequence, description: `Flanks CRISPR cut at ${guide.cut + 1}`, extra: {} });
      });
      state.panel = 'Primers';
      render();
    } catch (e) { toast(e.message, 'error'); }
  },
  async runFold({ start, end, maxLength }) {
    const doc = current();
    if (!doc) return;
    const seq = start != null && end != null ? subseq(doc.seq, start, end, doc.circular) : doc.seq;
    state.fold = { docId: doc.id, busy: true, maxLength };
    render();
    try {
      if (!seq.length) throw new Error('Nothing to fold');
      const { result, worker } = await offload(
        ({ sequence, limit }) => fold(sequence, { maxLength: limit }),
        { sequence: seq, limit: maxLength },
        { workerUrl: new URL('../workers/compute-worker.js', import.meta.url) },
      );
      if (result?.__error) throw new Error(result.__error);
      state.fold = { docId: doc.id, result, worker, maxLength, busy: false };
    } catch (e) {
      state.fold = { docId: doc.id, error: e.message, busy: false, maxLength };
    }
    render();
  },
  exportStructure() {
    const result = state.fold?.result;
    if (!result) return;
    const svg = structureToSvg(result.sequence, result.structure, { width: 900, perRow: 60 });
    const blob = new Blob([svg], { type: 'image/svg+xml' });
    const url = URL.createObjectURL(blob);
    const a = h('a', { href: url, download: 'structure.svg' });
    document.body.append(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 2000);
  },
  save(kind) { save(kind); },
  saveProject,
  openProject,
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
    box.append(h('div', { class: 'tab', role: 'tab', 'aria-selected': String(i === state.active), onclick: () => { state.active = i; state.selection = null; state.highlight = null; state.selectedFeatures = new Set(); state.traceBase = null; render(); el('panel').scrollTop = 0; } },
      h('span', { class: 'name' }, doc.name),
      dirty ? h('span', { class: 'dirty', title: 'edited' }, '•') : null,
      h('button', { class: 'close', title: 'Close', onclick: (e) => { e.stopPropagation(); closeDoc(i); } }, '×')));
  });
}

function renderViewBar(doc) {
  const opts = clear(el('viewOpts'));
  const ops = clear(el('editOps'));
  if (!doc) return;
  const toggles = [['showComplement', 'Both strands'], ['showTranslation', 'Translation'], ['showSixFrame', 'Six frames']];
  if (alignments.get(doc.id)?.length) toggles.push(['showReads', `Reads (${alignments.get(doc.id).length})`]);
  for (const [key, label] of toggles) {
    opts.append(h('label', { class: 'sub', style: { display: 'flex', gap: '4px', alignItems: 'center' } },
      h('input', { type: 'checkbox', checked: state[key], onchange: (e) => { state[key] = e.target.checked; render(); } }), label));
  }
  const sel = state.selection;
  ops.append(h('button', { class: 'btn', disabled: !sel, onclick: deleteSelection }, 'Delete'));
  ops.append(h('button', { class: 'btn', disabled: !sel, onclick: replaceSelection }, 'Replace…'));
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
  items.push({ label: 'Copy fragment with features', disabled: !sel, onClick: () => api.copyFragment() });
  items.push({ label: 'Paste fragment here', disabled: !state.clipboard, onClick: () => api.pasteFragment() });
  items.push('-');
  items.push({ label: 'Make a forward primer', disabled: !picked, onClick: () => addPrimerFromSelection(1) });
  items.push({ label: 'Make a reverse primer', disabled: !picked, onClick: () => addPrimerFromSelection(-1) });
  items.push({ label: 'Design a primer pair…', disabled: !picked, onClick: () => { state.panel = 'Primers'; render(); } });
  items.push('-');
  items.push({ label: 'Add a feature here', disabled: !sel, onClick: addFeatureFromSelection });
  items.push({ label: 'Delete selected bases', danger: true, disabled: !sel, onClick: deleteSelection });
  items.push({ label: 'Replace selected bases…', disabled: !sel, onClick: replaceSelection });
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
  activeTraceView = null;
  if (!doc) { view.append(welcome()); return; }
  if (!doc.seq.length && state.view !== 'trace') { view.append(emptyCard(doc)); return; }
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
      onSelect: (i) => { state.trace = i; state.traceBase = null; render(); },
      onZoom: (px) => { state.pxPerBase = px; render(); },
      highlightBase: state.traceBase,
      onPickBase: (k, trace) => {
        state.traceBase = k;
        const doc2 = current();
        if (doc2) {
          if (doc2.seq.length * trace.bases.length > 2_000_000) {
            toast('This trace is too large to map onto the reference automatically', 'error');
            return;
          }
          const key = `${doc2.id}:${trace.name ?? ''}:${trace.bases.length}`;
          let cached = traceMaps.get(doc2.id);
          if (!cached || cached.key !== key) {
            try {
              cached = { key, helper: traceBaseToReference(doc2.seq, trace, k) };
              traceMaps.set(doc2.id, cached);
            } catch (e) {
              toast(`Could not map the trace: ${e.message}`, 'error');
              return;
            }
          }
          const position = cached.helper.map[k];
          if (position != null) state.selection = { start: position, end: position + 1 };
          toast(`Trace base ${k + 1}: ${trace.bases[k]}${position != null ? ` → reference ${position + 1}` : ''}${trace.confidence ? ` · quality ${trace.confidence[k] ?? '?'}` : ''}`);
        }
        activeTraceView?.setHighlight?.(k);
        renderStatus(current());
        renderPanels(current());
      },
      onLocate: (k, trace) => {
        const doc2 = current();
        if (!doc2) return;
        if (doc2.seq.length * trace.bases.length > 2_000_000) { toast('This trace is too large to map onto the reference automatically', 'error'); return; }
        let position;
        try { position = traceBaseToReference(doc2.seq, trace, k).map[k]; } catch (e) { toast(`Could not map the trace: ${e.message}`, 'error'); return; }
        if (position == null) { toast('This trace base has no reference position', 'error'); return; }
        state.selection = { start: position, end: position + 1 };
        state.view = 'seq';
        render();
        toast(`Showing reference base ${position + 1}`);
      },
    });
    activeTraceView = tv;
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
        if (!info?.live) { syncHighlight(); render(); }
      },
      showTranslation: state.showTranslation, showSixFrame: state.showSixFrame, showComplement: state.showComplement,
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
      onclick: () => { state.panel = name; render(); el('panel').scrollTop = 0; },
    }, name));
  }
  const panel = clear(el('panel'));
  if (!doc) { panel.append(h('p', { class: 'empty' }, 'Open a file to see its details.')); return; }
  panel.append(renderPanel(state.panel, {
    doc, state, api,
    others: state.docs.filter((d) => d.id !== doc.id),
    docs: state.docs,
  }));
}

// A new document has no bases, and nothing to select, so give it a way in.
function emptyCard(doc) {
  const box = h('textarea', { id: 'pasteSeq', rows: 6, placeholder: 'ACGT… or a FASTA record', class: 'mono', style: { width: '100%' } });
  const add = () => {
    const text = box.value.replace(/^>.*$/gm, '');
    const seq = cleanSequence(text);
    if (!seq) { toast('No bases found in what was pasted', 'error'); return; }
    edit('Add sequence', (d) => { replaceRange(d, 0, 0, seq); });
  };
  return h('div', { class: 'welcome' },
    h('h1', {}, `${doc.name} is empty`),
    h('p', {}, 'Paste a sequence to start, or open a file. Spaces, digits and line breaks are ignored.'),
    h('label', { class: 'field' }, h('span', {}, 'Sequence'), box),
    h('div', { class: 'group' }, h('button', { class: 'btn primary', onclick: add }, 'Add sequence')));
}

function welcome() {
  if (state.restorable) return restoreCard(state.restorable);
  const drop = h('div', { class: 'drop' },
    h('p', { style: { margin: '0 0 12px' } }, 'Drop a sequence file here'),
    h('button', { class: 'btn primary', onclick: () => el('fileInput').click() }, 'Choose a file…'));
  const box = h('div', { class: 'welcome' },
    h('h1', {}, 'A plasmid map in the browser'),
    h('p', {}, 'Open a .dna, GenBank or FASTA file to see its map, sequence, features and cut sites. Files stay on this computer.'),
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

// Rebuilding the panels replaces every input in them. Remember which one the
// reader was typing in, and where the caret was, so a keystroke in the filter
// box does not end with the focus lost.
function focusSnapshot() {
  const a = document.activeElement;
  if (!a?.id || !(el('panel').contains(a) || el('viewOpts').contains(a) || el('view').contains(a))) return null;
  const range = typeof a.selectionStart === 'number' ? [a.selectionStart, a.selectionEnd, a.selectionDirection] : null;
  return { id: a.id, range };
}

function restoreFocus(snap) {
  if (!snap) return;
  const n = document.getElementById(snap.id);
  if (!n || n === document.activeElement) return;
  n.focus({ preventScroll: true });
  try { if (snap.range) n.setSelectionRange(...snap.range); } catch { /* inputs that have no caret, such as a range slider */ }
}

function render() {
  const doc = current();
  // drop references to features that were deleted or never existed in this document
  if (state.highlight && !doc?.features.some((f) => f.id === state.highlight)) state.highlight = null;
  if (state.selectedFeatures.size) {
    const ids = new Set(doc?.features.map((f) => f.id));
    const kept = [...state.selectedFeatures].filter((id) => ids.has(id));
    if (kept.length !== state.selectedFeatures.size) state.selectedFeatures = new Set(kept);
  }
  const snap = focusSnapshot();
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
  el('undoBtn').disabled = !(undo.get(doc?.id)?.canUndo?.());
  el('redoBtn').disabled = !(undo.get(doc?.id)?.canRedo?.());
  el('saveAs').disabled = !doc;
  document.title = doc ? `${doc.name} · BublaGene` : 'BublaGene';
  restoreFocus(snap);
  rememberSoon();
}

// ---- events -------------------------------------------------------------------

el('openBtn').addEventListener('click', async () => {
  if (supportsFileSystemAccess()) {
    try {
      const files = await openDocument({ multiple: true });
      if (files?.length) await openFiles(files);
      return;
    } catch (e) {
      if (e?.name !== 'AbortError') toast(`Could not open: ${e.message}`, 'error');
      return;
    }
  }
  el('fileInput').click();
});
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
loadSession().then(async (saved) => {
  if (saved?.docs?.length && !state.docs.length) {
    state.restorable = saved;
    render();
    return;
  }
  const backup = await readAutosave('bublagene.backup');
  if (!backup || state.docs.length) return;
  try {
    const parsed = JSON.parse(backup);
    if (parsed?.docs?.length) {
      state.restorable = { ...parsed, restoredFromBackup: true };
      render();
    }
  } catch { /* ignore a corrupt backup */ }
});
