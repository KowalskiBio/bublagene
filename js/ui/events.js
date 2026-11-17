// @ts-nocheck — legacy DOM wiring; the engines it calls are type-checked.
// Buttons, keys, drag-and-drop, clipboard and print events.

import { newDocument } from '../model.js';
import { cleanSequence } from '../seq.js';
import { openDocument, supportsFileSystemAccess } from '../storage.js';
import { api } from './api/index.js';
import { copy, selectedSeq } from './clipboard.js';
import { openContextMenu } from './contextmenu.js';
import { current, el, render, state } from './core.js';
import { openDoc, undoOnce } from './documents.js';
import { fmt, toast } from './dom.js';
import { printSheet } from './export.js';
import { openFiles, save } from './files.js';
import { closeMenu } from './menu.js';
import { remember } from './persist.js';
import { BASE_KEY, deleteAtCaret, insertAtCaret, moveCaret } from './typing.js';
import { overlaysFor, shownSites } from './viewdata.js';
import { renderView } from './viewpane.js';

function buildPrintSheet() {
  const doc = current();
  const root = document.getElementById('printRoot');
  if (!root || !doc) return;
  root.replaceChildren(printSheet(doc, { sites: shownSites(doc), style: state.mapStyle, overlays: overlaysFor(doc) }));
}

// Attach every listener the page needs. Called once, after the first render.
export function bindEvents() {
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
    else if (!typing && e.key === 'Escape') { closeMenu(); state.selection = null; state.caret = null; state.highlight = null; render(); }
    else if (!typing && !mod && !e.altKey && current() && state.view === 'seq') {
      // typing straight into the sequence view
      let handled = false;
      if (BASE_KEY.test(e.key)) handled = insertAtCaret(e.key.toUpperCase());
      else if (e.key === 'Backspace') handled = deleteAtCaret(-1);
      else if (e.key === 'Delete') handled = deleteAtCaret(1);
      else if (['ArrowLeft', 'ArrowRight', 'Home', 'End'].includes(e.key)) handled = moveCaret(e.key, e.shiftKey);
      if (handled) e.preventDefault();
    }
  });

  // pasting bases at the caret, or over the selection
  document.addEventListener('paste', (e) => {
    if (/^(INPUT|TEXTAREA|SELECT)$/.test(e.target.tagName) || e.target.isContentEditable) return;
    if (!current() || state.view !== 'seq') return;
    const text = e.clipboardData?.getData('text/plain') ?? '';
    const seq = cleanSequence(text.replace(/^>.*$/gm, ''));
    if (!seq) return;
    e.preventDefault();
    if (!insertAtCaret(seq, 'Paste bases')) toast('Click in the sequence first, or select the bases to replace', 'error');
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

  window.addEventListener('resize', () => { if (state.view === 'seq' || state.view === 'split') renderView(current()); });

  // The printed page is a sheet built for paper, complete however little of the
  // sequence is on screen; it exists only while printing.
  window.addEventListener('beforeprint', buildPrintSheet);
  window.addEventListener('afterprint', () => { document.getElementById('printRoot')?.replaceChildren(); });

  // write the session out before the tab goes away, rather than waiting out the debounce
  addEventListener('pagehide', remember);
  document.addEventListener('visibilitychange', () => { if (document.visibilityState === 'hidden') remember(); });
}
