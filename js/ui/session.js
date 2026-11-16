// Remembers the open documents between visits.
//
// IndexedDB rather than localStorage: a .dna file's untouched segments are kept
// as raw bytes so they can be written back, and IndexedDB stores typed arrays
// as they are instead of turning them into JSON. Every call is wrapped, because
// storage can be unavailable (a private window, blocked site data) and that
// must never stop the app from working.

const DB_NAME = 'bublagene';
const STORE = 'session';
const KEY = 'current';
const VERSION = 1;

let dbPromise = null;

function openDb() {
  if (dbPromise) return dbPromise;
  dbPromise = new Promise((resolve, reject) => {
    if (!globalThis.indexedDB) { reject(new Error('no IndexedDB')); return; }
    const req = indexedDB.open(DB_NAME, VERSION);
    req.onupgradeneeded = () => {
      if (!req.result.objectStoreNames.contains(STORE)) req.result.createObjectStore(STORE);
    };
    req.onsuccess = () => resolve(req.result);
    req.onerror = () => reject(req.error ?? new Error('IndexedDB refused'));
    req.onblocked = () => reject(new Error('IndexedDB blocked'));
  }).catch((e) => { dbPromise = null; throw e; });
  return dbPromise;
}

function run(mode, fn) {
  return openDb().then((db) => new Promise((resolve, reject) => {
    const tx = db.transaction(STORE, mode);
    const req = fn(tx.objectStore(STORE));
    tx.oncomplete = () => resolve(req?.result);
    tx.onerror = () => reject(tx.error);
    tx.onabort = () => reject(tx.error);
  }));
}

export async function saveSession(payload) {
  try {
    await run('readwrite', (store) => store.put(payload, KEY));
    return true;
  } catch {
    return false;
  }
}

export async function loadSession() {
  try {
    return await run('readonly', (store) => store.get(KEY));
  } catch {
    return null;
  }
}

export async function clearSession() {
  try {
    await run('readwrite', (store) => store.delete(KEY));
    return true;
  } catch {
    return false;
  }
}

export function timeAgo(iso) {
  const then = new Date(iso).getTime();
  if (!Number.isFinite(then)) return 'earlier';
  const mins = Math.round((Date.now() - then) / 60000);
  if (mins < 1) return 'a moment ago';
  if (mins === 1) return 'a minute ago';
  if (mins < 60) return `${mins} minutes ago`;
  const hours = Math.round(mins / 60);
  if (hours === 1) return 'an hour ago';
  if (hours < 24) return `${hours} hours ago`;
  const days = Math.round(hours / 24);
  return days === 1 ? 'yesterday' : `${days} days ago`;
}
