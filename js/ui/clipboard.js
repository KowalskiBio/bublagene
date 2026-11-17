// @ts-nocheck — legacy DOM wiring; the engines it calls are type-checked.
// Copying bases, and turning the selection into primers and features.

import { colorFor, newId } from '../model.js';
import { revcomp, subseq } from '../seq.js';
import { current, state } from './core.js';
import { edit } from './documents.js';
import { fmt, toast } from './dom.js';

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

export function copy(text, what = 'bases') {
  if (!text) return;
  const done = `Copied ${fmt(text.length)} ${what}`;
  if (copyViaTextarea(text)) { toast(done); return; }
  navigator.clipboard?.writeText(text).then(
    () => toast(done),
    () => toast('Could not reach the clipboard', 'error'),
  ) ?? toast('Could not reach the clipboard', 'error');
}

export function selectedSeq() {
  const doc = current();
  if (!doc || !state.selection) return '';
  return subseq(doc.seq, state.selection.start, state.selection.end, doc.circular);
}

export function addPrimerFromSelection(dir) {
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

export function addFeatureFromSelection() {
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
