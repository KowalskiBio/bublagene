// @ts-nocheck — legacy DOM wiring; the engines it calls are type-checked.
// Opening, closing and editing documents, with undo.

import { ENZYME_BY_NAME, findSites } from '../enzymes.js';
import { readAlignments } from '../formats/bam.js';
import { readTraces } from '../formats/ztr.js';
import { placeRead } from '../readalign.js';
import { readTrace, traceToDocument } from '../traces.js';
import { addHistory, newId, cloneDoc, featureSpan, replaceRange } from '../model.js';
import { cleanSequence } from '../seq.js';
import { History } from '../undo.js';
import { alignments, current, isDirty, render, state, traces, undo } from './core.js';
import { toast } from './dom.js';

export async function loadAlignment(doc) {
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

// Add an aligned sequence to a reference document, with its chromatogram when
// it has one, and show it in the Alignment view. `read` carries `map` (the
// reference position of each base of the sequence).
export function attachRead(ref, read, trace = null) {
  const { map, ...row } = read;
  if (trace) {
    trace.reverse = read.reverse;
    trace.refMap = map;
    row.traceId = trace.id;
    traces.set(ref.id, [...(traces.get(ref.id) ?? []).filter((t) => t.id !== trace.id), trace]);
  }
  const reads = [...(alignments.get(ref.id) ?? []), { ...row, pos: row.start }]
    .sort((a, b) => a.start - b.start || a.label.localeCompare(b.label));
  alignments.set(ref.id, reads);
  state.showReads = true;
  state.alignRow = reads.findIndex((r) => r.label === row.label && r.start === row.start);
  state.alignScroll = null;
  state.view = 'align';
  state.active = state.docs.indexOf(ref);
  render();
}

// Open a chromatogram file. With a reference open (a .gb or .dna the user
// opened first) the read is aligned to it and added as a row under the
// sequence, with its chromatogram behind the Trace view. Otherwise, or if the
// read does not match, it opens as a document of its own.
export async function openTraceFile(buffer, fileName) {
  const trace = { id: newId(), ...(await readTrace(buffer, fileName)) };
  if (!trace.name) trace.name = fileName.replace(/\.[^.]+$/, '');
  const ref = current();
  if (ref && ref.seq.length && ref.source?.format !== 'trace') {
    const read = placeRead(ref.seq, trace, { label: trace.name });
    if (read) {
      attachRead(ref, read, trace);
      toast(`${trace.name}: aligned to ${ref.name} at ${read.start + 1}–${read.end}, ${read.identity.toFixed(1)}% identity${read.reverse ? ', reverse strand' : ''}`);
      return;
    }
    toast(`${trace.name} does not match ${ref.name}; opened on its own`, 'error');
  }
  const doc = traceToDocument(trace);
  openDoc(doc);
  traces.set(doc.id, [trace]);
  state.view = 'trace';
  state.trace = 0;
  render();
}

export function openDoc(doc) {
  state.restorable = null; // opening something is an answer to the restore offer
  state.docs.push(doc);
  state.active = state.docs.length - 1;
  undo.set(doc.id, new History());
  state.selection = null;
  state.caret = null;
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

export function closeDoc(i) {
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
    state.caret = null;
    state.highlight = null;
    state.selectedFeatures = new Set();
    state.traceBase = null;
  }
  render();
}

// Callers that adjust the selection or caret afterwards pass { redraw: false } and
// render once themselves, so a keystroke does not redraw the page twice.
export function edit(label, mutate, { redraw = true } = {}) {
  const doc = current();
  if (!doc) return;
  const stack = undo.get(doc.id);
  stack.record(cloneDoc(doc), label, { coalesceKey: label });
  mutate(doc);
  state.alignResult = null;
  state.fold = null;
  state.crispr = null;
  state.detect = null;
  state.mutate = { ...state.mutate, result: null };
  // a selection that no longer fits the edited sequence would point past its end
  if (state.selection && (state.selection.end > doc.seq.length || state.selection.start > doc.seq.length)) state.selection = null;
  if (state.caret !== null && state.caret > doc.seq.length) state.caret = doc.seq.length;
  addHistory(doc, 'edit', label);
  if (redraw) render();
}

// Replace the selection with typed text. Cancelling the prompt must leave the
// sequence alone; an empty answer is treated the same way (Delete is for that).
export function replaceSelection() {
  const sel = state.selection;
  if (!sel) return;
  const answer = window.prompt('Replace the selection with:');
  if (answer === null) return;
  const seq = cleanSequence(answer);
  if (!seq) return;
  edit('Replace selection', (d) => { replaceRange(d, sel.start, sel.end, seq); }, { redraw: false });
  state.selection = { start: sel.start, end: sel.start + seq.length };
  render();
}

export function deleteSelection() {
  const sel = state.selection;
  if (!sel) return;
  edit('Delete selection', (d) => { replaceRange(d, sel.start, sel.end, ''); }, { redraw: false });
  state.selection = null;
  render();
}

// A highlighted feature is only "current" while the selection is that feature;
// once the reader selects something else, or the feature is gone, drop it.
export function syncHighlight() {
  if (!state.highlight) return;
  const f = current()?.features.find((x) => x.id === state.highlight);
  const sel = state.selection;
  const span = f && featureSpan(f);
  if (!f || !sel || sel.start !== span.start || sel.end !== span.end) state.highlight = null;
}

export function undoOnce(redo = false) {
  const doc = current();
  if (!doc) return;
  const stack = undo.get(doc.id);
  const prev = redo ? stack.redo(cloneDoc(doc)) : stack.undo(cloneDoc(doc));
  if (!prev) return;
  state.docs[state.active] = prev;
  undo.set(prev.id, stack);
  state.selection = null;
  state.caret = null;
  state.selectedFeatures = new Set();
  render();
}
