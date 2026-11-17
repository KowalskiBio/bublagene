// Sequence-keyed caches and stale-record invalidation.
//
// Caches and records store the digest of the sequence they describe. Any edit
// changes the digest, so a lookup either misses or marks the old record stale.

export function seqDigest(seq) {
  let h1 = 0x811c9dc5;
  let h2 = 0x9e3779b9;
  const s = String(seq ?? '');
  for (let i = 0; i < s.length; i++) {
    const c = s.charCodeAt(i) & 0xff;
    h1 = Math.imul(h1 ^ c, 0x01000193) >>> 0;
    h2 = Math.imul(h2 ^ (c + i), 0x85ebca6b) >>> 0;
  }
  return `${h1.toString(16).padStart(8, '0')}${h2.toString(16).padStart(8, '0')}`;
}

export class SequenceCache {
  constructor({ limit = 64 } = {}) {
    this.limit = limit;
    this.entries = new Map();
  }

  get(name, seq, compute) {
    const digest = seqDigest(seq);
    const hit = this.entries.get(name);
    if (hit && hit.digest === digest) return hit.value;
    const value = compute(seq);
    this.entries.set(name, { digest, value });
    if (this.entries.size > this.limit) this.entries.delete(this.entries.keys().next().value);
    return value;
  }

  has(name, seq) {
    const hit = this.entries.get(name);
    return Boolean(hit && hit.digest === seqDigest(seq));
  }

  invalidate(name) { this.entries.delete(name); }
  clear() { this.entries.clear(); }
}

// Alignment records keep the digest of the sequence they were computed against.
export function markAlignmentStale(doc) {
  if (!doc?.alignment) return doc;
  const digest = seqDigest(doc.seq);
  doc.alignment.stale = Boolean(doc.alignment.digest && doc.alignment.digest !== digest);
  doc.alignment.currentDigest = digest;
  return doc;
}

export function alignmentIsCurrent(doc) {
  if (!doc?.alignment) return true;
  return doc.alignment.digest === seqDigest(doc.seq);
}

export function stampAlignment(doc, record) {
  if (!doc || !record) return record;
  doc.alignment = { ...record, digest: seqDigest(doc.seq), currentDigest: seqDigest(doc.seq), stale: false };
  return doc.alignment;
}

// Named indexes (ORFs, enzyme sites, ...) recomputed only when the digest moves.
export function getIndex(doc, name, compute, { cache = doc.__cache ?? (doc.__cache = new SequenceCache()) } = {}) {
  return cache.get(name, doc.seq, compute);
}

/**
 * Remember the latest answer of a whole-sequence scan. The first argument is the
 * sequence; `keyOf` turns the other arguments into a string. One entry is kept per
 * key, so memory stays bounded and an edit simply replaces the stale answer.
 * Callers must treat the returned value as read-only.
 * @template {(seq: string, ...rest: any[]) => any} F
 * @param {F} fn
 * @param {(...rest: any[]) => string} keyOf
 * @param {number} [limit]
 * @returns {F}
 */
export function memoLatest(fn, keyOf, limit = 1024) {
  /** @type {Map<string, {seq: string, value: any}>} */
  const store = new Map();
  return /** @type {F} */ ((seq, ...rest) => {
    const key = keyOf(...rest);
    const hit = store.get(key);
    if (hit && hit.seq === seq) return hit.value;
    const value = fn(seq, ...rest);
    store.delete(key);
    store.set(key, { seq, value });
    if (store.size > limit) store.delete(store.keys().next().value);
    return value;
  });
}
