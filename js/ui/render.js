// @ts-nocheck — legacy DOM wiring; the engines it calls are type-checked.
// The one function that redraws the page. Importing this module registers it.

import { renderPanels, renderStatus, renderTabs, renderViewBar } from './chrome.js';
import { alignments, current, el, setRenderer, state, traces, undo } from './core.js';
import { rememberSoon } from './persist.js';
import { renderView } from './viewpane.js';

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

function renderApp() {
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
  const hasAligned = !!alignments.get(doc?.id)?.length;
  for (const b of el('viewSeg').querySelectorAll('button')) {
    b.setAttribute('aria-pressed', String(b.dataset.view === state.view));
    if (b.dataset.view === 'trace') b.hidden = !hasTraces;
    if (b.dataset.view === 'align') b.hidden = !hasAligned;
  }
  if (state.view === 'trace' && !hasTraces) state.view = 'map';
  if (state.view === 'align' && !hasAligned) state.view = 'map';
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

setRenderer(renderApp);
