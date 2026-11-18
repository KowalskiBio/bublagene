// @ts-nocheck — legacy DOM wiring; the engines it calls are type-checked.
// The main view: map, sequence, both together, or the chromatogram.

import { traceBaseToReference } from '../traces.js';
import { api } from './api/index.js';
import { emptyCard, welcome } from './cards.js';
import { renderPanels, renderStatus } from './chrome.js';
import { current, el, render, state, traces } from './core.js';
import { syncHighlight } from './documents.js';
import { clear, h, toast } from './dom.js';
import { circularMap, linearMap } from './map.js';
import { isOverlay, isPrimer } from './overlays.js';
import { sequenceView } from './seqview.js';
import { traceView } from './traceview.js';
import { overlaysFor, readsFor, shownSites } from './viewdata.js';

// the cache of trace-to-reference mappings, per document
const traceMaps = new Map();
let activeTraceView = null;

export function renderView(doc) {
  // keep the reader where they were: rebuilding the view resets its scroll
  const wasScrolled = state.view === 'seq' ? el('view').scrollTop : 0;
  const view = clear(el('view'));
  activeTraceView = null;
  if (!doc) { view.append(welcome()); return; }
  if (!doc.seq.length && state.view !== 'trace') { view.append(emptyCard(doc)); return; }
  const sites = shownSites(doc);
  const opts = {
    selection: state.selection, sites, highlight: state.highlight,
    overlays: overlaysFor(doc),
    // a reading frame is not a feature of the document: picking one just selects it
    onPick: (f) => {
      if (!f) return;
      if (!isOverlay(f)) { api.focusFeature(f); return; }
      if (isPrimer(f)) state.panel = 'Primers';
      api.setSelection({ start: f.segments[0].start, end: f.segments[0].end });
    },
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
        if (doc2 && trace.refMap) {
          // a read placed on this reference already knows where each base sits
          const position = trace.refMap[k];
          if (position != null) state.selection = { start: position, end: position + 1 };
          toast(`Trace base ${k + 1}: ${trace.bases[k]}${position != null ? ` → reference ${position + 1}` : ''}${trace.confidence ? ` · quality ${trace.confidence[k] ?? '?'}` : ''}`);
        } else if (doc2) {
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
        if (trace.refMap) {
          const at = trace.refMap[k];
          if (at == null) { toast('This trace base has no reference position', 'error'); return; }
          state.selection = { start: at, end: at + 1 };
          state.view = 'seq';
          render();
          toast(`Showing reference base ${at + 1}`);
          return;
        }
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
  const mapNode = () => ((state.mapStyle === 'auto' ? doc.circular : state.mapStyle === 'circular') ? circularMap : linearMap)(doc, opts);
  const mapBox = (extra = {}) => {
    const wrap = h('div', { class: 'mapwrap', ...extra });
    wrap.append(mapNode());
    wrap.addEventListener('click', (e) => { if (e.target === wrap) { state.highlight = null; state.selection = null; state.caret = null; render(); } });
    return wrap;
  };
  if (state.view === 'map') {
    view.append(mapBox());
    return;
  }
  const buildSequence = () => sequenceView(doc, {
    ...opts,
    // While the pointer is down the view repaints itself; rebuilding it here
    // would drop the drag. Only once the drag ends do the panels catch up.
    onSelect: (sel, info) => {
      state.selection = sel;
      state.selAnchor = null;
      state.selHead = null;
      if (sel) state.caret = null;
      else if (info && 'caret' in info && info.caret !== null) state.caret = info.caret;
      renderStatus(doc);
      if (!info?.live) { syncHighlight(); render(); }
    },
    showTranslation: state.showTranslation, showSixFrame: state.showSixFrame, showComplement: state.showComplement,
    reads: readsFor(doc), caret: state.selection ? null : state.caret,
  });
  if (state.view === 'split') {
    // the map on top, the sequence beneath it in its own scrolling box
    const box = h('div', { class: 'split' });
    const scroller = h('div', { class: 'splitseq' });
    const sv = buildSequence();
    scroller.append(sv);
    box.append(mapBox({ class: 'mapwrap splitmap' }), scroller);
    view.append(box);
    scroller.onscroll = () => { state.splitScroll = scroller.scrollTop; sv.onScroll?.(); };
    if (state.splitScroll) requestAnimationFrame(() => { scroller.scrollTop = state.splitScroll; sv.onScroll?.(); });
    return;
  }
  const sv = buildSequence();
  view.append(sv);
  view.onscroll = () => sv.onScroll?.();
  if (wasScrolled) {
    requestAnimationFrame(() => { view.scrollTop = wasScrolled; sv.onScroll?.(); });
  }
}
