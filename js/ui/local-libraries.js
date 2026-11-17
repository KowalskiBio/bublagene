// @ts-nocheck — legacy DOM wiring; the engines it calls are type-checked.
// Named enzyme sets and detection elements the reader saved. They live in this browser only.

import { normalizeLibrary } from '../detect.js';
import { ENZYME_BY_NAME } from '../enzymes.js';

// Named enzyme sets the reader has saved. They live in this browser only.
const SETS_KEY = 'bublagene.enzymeSets';
export function loadEnzymeSets() {
  try {
    const raw = JSON.parse(localStorage.getItem(SETS_KEY) ?? '{}');
    const out = {};
    for (const [name, list] of Object.entries(raw ?? {})) {
      if (Array.isArray(list)) out[name] = list.filter((n) => ENZYME_BY_NAME[n]);
    }
    return out;
  } catch { return {}; }
}
export function storeEnzymeSets(sets) {
  try { localStorage.setItem(SETS_KEY, JSON.stringify(sets)); return true; } catch { return false; }
}

// Elements the reader saved for detection; kept in this browser only.
const LIB_KEY = 'bublagene.featureLibrary';
export function loadUserLibrary() {
  try { return normalizeLibrary(JSON.parse(localStorage.getItem(LIB_KEY) ?? '[]')); } catch { return []; }
}
export function storeUserLibrary(list) {
  try { localStorage.setItem(LIB_KEY, JSON.stringify(list)); return true; } catch { return false; }
}
