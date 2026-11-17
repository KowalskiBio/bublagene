// @ts-nocheck — legacy DOM wiring; the engines it calls are type-checked.
// Remembering the session between visits, and offering to reopen it.

import { COMMON_SET } from '../enzymes.js';
import { cloneDoc } from '../model.js';
import { readAutosave, writeAutosave } from '../storage.js';
import { History } from '../undo.js';
import { isDirty, render, state, undo } from './core.js';
import { loadAlignment } from './documents.js';
import { toast } from './dom.js';
import { clearSession, loadSession, saveSession } from './session.js';

let saveTimer = null;
export function rememberSoon() {
  if (state.restorable) return; // nothing of the user's own is open yet
  clearTimeout(saveTimer);
  saveTimer = setTimeout(remember, 900);
}

export function remember() {
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
    showOrfs: state.showOrfs,
    mapStyle: state.mapStyle,
    tm: state.tm,
    showPrimers: state.showPrimers,
    orfMin: state.orfMin,
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

export function restoreSession(saved) {
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
  state.showOrfs = saved.showOrfs ?? false;
  state.mapStyle = saved.mapStyle ?? 'auto';
  state.tm = saved.tm ?? state.tm;
  state.showPrimers = saved.showPrimers ?? true;
  state.orfMin = saved.orfMin ?? 100;
  state.showComplement = saved.showComplement ?? true;
  state.featureFilter = saved.featureFilter ?? '';
  state.restorable = null;
  render();
  for (const d of state.docs) loadAlignment(d);
  toast(`Reopened ${state.docs.length === 1 ? state.docs[0].name : `${state.docs.length} files`}`);
}

export function startFresh() {
  state.restorable = null;
  clearSession();
  render();
}

// offer the previous session once storage answers
export function offerPreviousSession() {
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
}
