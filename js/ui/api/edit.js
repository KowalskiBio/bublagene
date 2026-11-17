// @ts-nocheck — legacy DOM wiring; the engines it calls are type-checked.
// The panels' view of the app: Selections, features, the clipboard, fragments and the undo timeline.

import { copyFeatures, extractFragment, pasteFeatures, pasteFragment, reorderFeature } from '../../fragment.js';
import { cloneDoc, colorFor, featureSequence, featureSpan, newId } from '../../model.js';
import { addPrimerFromSelection, copy, selectedSeq } from '../clipboard.js';
import { current, el, render, state, undo } from '../core.js';
import { edit, openDoc, syncHighlight } from '../documents.js';
import { fmt, toast } from '../dom.js';

export function editApi(api) {
  return {
    edit,
    openDoc,
    setSelection(sel) { state.selection = sel; state.selAnchor = null; state.selHead = null; if (sel) state.caret = null; syncHighlight(); render(); },
    setFeatureFilter(text) { state.featureFilter = text; render(); },
    focusFeature(f) {
      state.highlight = state.highlight === f.id ? null : f.id;
      if (state.highlight) {
        state.selection = { start: featureSpan(f).start, end: featureSpan(f).end };
        state.panel = 'Features'; // so the details appear in the sidebar
      }
      render();
    },
    primerFromSelection(dir) { addPrimerFromSelection(dir); },
    selectedSequence() { return selectedSeq(); },
    copy(text, what) { copy(text, what); },
    addOrfAsCds(o) {
      const id = newId();
      state.panel = 'Features';
      state.highlight = id;
      state.selection = { start: o.start, end: o.end };
      edit(`Add CDS ${o.start + 1}–${o.end || current().seq.length}`, (d) => {
        d.features.push({
          id, name: `ORF ${o.codons} aa`, type: 'CDS', strand: o.strand, color: colorFor('CDS'),
          segments: [{ start: o.start, end: o.end }], qualifiers: [], extra: {},
        });
      });
    },
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
    copyFeatureAnnotations(ids = null) {
      const doc = current();
      if (!doc) return false;
      const pick = ids ?? (state.selectedFeatures.size ? [...state.selectedFeatures] : state.highlight ? [state.highlight] : []);
      const clip = copyFeatures(doc, pick);
      if (!clip) { toast('Tick the features to copy first', 'error'); return false; }
      state.featureClip = clip;
      toast(`Copied ${clip.features.length} feature${clip.features.length === 1 ? '' : 's'}`);
      render();
      return true;
    },
    cutFeatureAnnotations() {
      const doc = current();
      if (!doc) return;
      const ids = state.selectedFeatures.size ? [...state.selectedFeatures] : state.highlight ? [state.highlight] : [];
      if (!api.copyFeatureAnnotations(ids)) return;
      const gone = new Set(ids);
      edit('Cut features', (d) => { d.features = d.features.filter((f) => !gone.has(f.id)); });
      state.selectedFeatures = new Set();
      render();
    },
    pasteFeatureAnnotations() {
      const doc = current();
      const clip = state.featureClip;
      if (!doc || !clip) return;
      const at = state.selection?.start ?? state.caret ?? clip.origin;
      let result;
      try {
        edit('Paste features', (d) => { result = pasteFeatures(d, at, clip); });
      } catch (e) { toast(e.message, 'error'); return; }
      toast(result.skipped
        ? `Pasted ${result.added}; ${result.skipped} would run past the end`
        : `Pasted ${result.added} feature${result.added === 1 ? '' : 's'}`);
    },
    copyFragment() {
      const doc = current();
      if (!doc || !state.selection) return;
      state.clipboard = extractFragment(doc, state.selection);
      if (state.clipboard) toast(`Copied ${fmt(state.clipboard.seq.length)} bp with ${state.clipboard.features.length} feature(s)`);
    },
    pasteFragment() {
      const doc = current();
      if (!doc || !state.clipboard) return;
      const at = state.selection?.start ?? state.caret ?? 0;
      edit('Paste fragment', (d) => pasteFragment(d, at, state.clipboard));
      toast(`Pasted ${fmt(state.clipboard.seq.length)} bp`);
    },
    focusDoc(id) {
      const i = state.docs.findIndex((d) => d.id === id);
      if (i < 0) return;
      state.active = i; state.selection = null; state.caret = null; state.highlight = null; state.selectedFeatures = new Set();
      render(); el('panel').scrollTop = 0;
    },
    goTo(docId, start, end) {
      const i = state.docs.findIndex((d) => d.id === docId);
      if (i < 0) return;
      if (i !== state.active) { state.active = i; state.highlight = null; state.selectedFeatures = new Set(); state.traceBase = null; }
      state.caret = null;
      state.selection = end > start ? { start, end } : null;
      if (state.view === 'map') state.view = 'seq';
      render();
    },
    // ---- history ----
    timeline() { const d = current(); return (d && undo.get(d.id)?.timeline?.()) ?? { labels: [], position: 0 }; },
    jumpHistory(position) {
      const doc = current();
      if (!doc) return;
      const stack = undo.get(doc.id);
      const prev = stack.jumpTo(position, cloneDoc(doc));
      state.docs[state.active] = prev;
      undo.set(prev.id, stack);
      state.selection = null; state.caret = null; state.highlight = null; state.selectedFeatures = new Set();
      render();
    },
  };
}
