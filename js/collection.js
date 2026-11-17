// Collections: files the reader chose to keep in this browser, grouped by name,
// and searchable without opening them. IndexedDB, separate from the session
// store, so clearing a session never touches a saved file. Every call can
// reject (private windows, blocked storage); callers show that and carry on.

import { cloneDoc } from './model.js';

const DB_NAME = 'bublagene-library';
const STORE = 'docs';

let dbPromise = null;

function openDb() {
  if (dbPromise) return dbPromise;
  dbPromise = new Promise((resolve, reject) => {
    if (!globalThis.indexedDB) { reject(new Error('This browser has no IndexedDB')); return; }
    const req = indexedDB.open(DB_NAME, 1);
    req.onupgradeneeded = () => {
      if (!req.result.objectStoreNames.contains(STORE)) req.result.createObjectStore(STORE, { keyPath: 'id' });
    };
    req.onsuccess = () => resolve(req.result);
    req.onerror = () => reject(req.error ?? new Error('IndexedDB refused'));
    req.onblocked = () => reject(new Error('IndexedDB blocked'));
  }).catch((e) => { dbPromise = null; throw e; });
  return dbPromise;
}

const run = (mode, fn) => openDb().then((db) => new Promise((resolve, reject) => {
  const tx = db.transaction(STORE, mode);
  const req = fn(tx.objectStore(STORE));
  tx.oncomplete = () => resolve(req?.result);
  tx.onerror = () => reject(tx.error);
  tx.onabort = () => reject(tx.error);
}));

// What a record needs to be listed and searched without loading the whole file.
export function recordFor(doc, collection = 'My files') {
  return {
    id: doc.id,
    collection: collection.trim() || 'My files',
    name: doc.name,
    savedAt: new Date().toISOString(),
    bp: doc.seq.length,
    circular: !!doc.circular,
    featureNames: doc.features.map((f) => f.name),
    primerNames: (doc.primers ?? []).map((p) => p.name),
    description: doc.notes?.description ?? '',
    doc: cloneDoc(doc),
  };
}

export const saveToCollection = (doc, collection) => run('readwrite', (s) => s.put(recordFor(doc, collection)));
export const removeFromCollection = (id) => run('readwrite', (s) => s.delete(id));
export const getFromCollection = (id) => run('readonly', (s) => s.get(id)).then((r) => r?.doc ?? null);

export async function listCollection() {
  const all = await run('readonly', (s) => s.getAll());
  return (all ?? []).map(({ doc, ...meta }) => meta).sort((a, b) => a.collection.localeCompare(b.collection) || a.name.localeCompare(b.name));
}

// Find files whose name, features, primers or description mention `query`, or
// whose sequence contains it (either strand). `loadSequences` reads whole files
// only when the query looks like DNA.
export async function searchCollection(query, { sequence = true } = {}) {
  const q = String(query ?? '').trim();
  if (!q) return [];
  const lower = q.toLowerCase();
  const all = await run('readonly', (s) => s.getAll());
  const isDna = /^[ACGTUNRYKMSWBDHV]{6,}$/i.test(q);
  const results = [];
  for (const r of all ?? []) {
    const reasons = [];
    if (r.name.toLowerCase().includes(lower)) reasons.push('name');
    const f = r.featureNames.filter((n) => n.toLowerCase().includes(lower));
    if (f.length) reasons.push(`feature ${[...new Set(f)].slice(0, 3).join(', ')}`);
    const p = r.primerNames.filter((n) => n.toLowerCase().includes(lower));
    if (p.length) reasons.push(`primer ${[...new Set(p)].slice(0, 3).join(', ')}`);
    if ((r.description ?? '').toLowerCase().includes(lower)) reasons.push('description');
    if (sequence && isDna && r.doc?.seq) {
      const up = q.toUpperCase().replace(/U/g, 'T');
      const rc = [...up].reverse().map((c) => ({ A: 'T', C: 'G', G: 'C', T: 'A' }[c] ?? c)).join('');
      const hay = r.circular ? r.doc.seq + r.doc.seq.slice(0, up.length - 1) : r.doc.seq;
      const fwd = hay.indexOf(up);
      const rev = hay.indexOf(rc);
      if (fwd >= 0) reasons.push(`sequence at ${fwd + 1}`);
      else if (rev >= 0) reasons.push(`sequence (reverse strand) at ${rev + 1}`);
    }
    if (reasons.length) { const { doc, ...meta } = r; results.push({ ...meta, reasons }); }
  }
  return results;
}
