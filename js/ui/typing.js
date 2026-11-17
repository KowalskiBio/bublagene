// @ts-nocheck — legacy DOM wiring; the engines it calls are type-checked.
// Typing straight into the sequence view: insert, delete and move the caret.

import { replaceRange } from '../model.js';
import { current, render, state } from './core.js';
import { deleteSelection, edit, syncHighlight } from './documents.js';

export const BASE_KEY = /^[acgtunrykmswbdhv]$/i;

// Put text where the reader is: over the selection if there is one, else at
// the caret. Returns false when there is nowhere to put it.
export function insertAtCaret(text, label = 'Type bases') {
  const doc = current();
  const sel = state.selection;
  if (!doc || !text) return false;
  const at = sel ? sel.start : state.caret;
  if (at === null || at === undefined) return false;
  const wraps = !!sel && sel.end < sel.start; // a wrapping selection rotates the origin to its start
  edit(label, (d) => { replaceRange(d, at, sel ? sel.end : at, text); }, { redraw: false });
  state.selection = null;
  state.caret = (wraps ? 0 : at) + text.length;
  render();
  return true;
}

// Backspace (dir -1) and Delete (dir +1)
export function deleteAtCaret(dir) {
  const doc = current();
  if (!doc) return false;
  if (state.selection) { deleteSelection(); return true; }
  const at = state.caret;
  if (at === null) return false;
  const from = dir < 0 ? at - 1 : at;
  if (from < 0 || from + 1 > doc.seq.length) return true; // at an end: nothing to remove
  edit('Delete bases', (d) => { replaceRange(d, from, from + 1, ''); }, { redraw: false });
  state.caret = from;
  render();
  return true;
}

// Arrow keys move the insertion point; with Shift they grow a selection from it.
export function moveCaret(key, shift) {
  const doc = current();
  if (!doc) return false;
  const L = doc.seq.length;
  const sel = state.selection;
  const step = key === 'ArrowLeft' ? -1 : key === 'ArrowRight' ? 1 : key === 'Home' ? -Infinity : Infinity;
  if (shift) {
    // the moving edge is the one away from where the selection began
    const anchor = state.selAnchor ?? (sel ? sel.start : state.caret);
    if (anchor === null || anchor === undefined) return false;
    const head = Math.max(0, Math.min(L, (state.selHead ?? (sel ? sel.end : state.caret)) + step));
    state.selAnchor = anchor;
    state.selHead = head;
    state.selection = head === anchor ? null : { start: Math.min(anchor, head), end: Math.max(anchor, head) };
    state.caret = state.selection ? null : head;
  } else {
    let at;
    if (sel) at = step < 0 ? sel.start : sel.end;
    else if (state.caret !== null) at = Math.max(0, Math.min(L, state.caret + step));
    else return false;
    state.selection = null;
    state.selAnchor = null;
    state.selHead = null;
    state.caret = at;
  }
  syncHighlight();
  render();
  return true;
}
