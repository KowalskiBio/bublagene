// Automatic feature detection: find known elements in a sequence.
//
// A library entry is { name, type, kind: 'dna' | 'protein', seq }. DNA entries
// are matched on both strands with a mismatch budget; protein entries are
// matched in all six frames. Wrap-around hits are found on circular sequences.

import { revcomp, translate, rangeLength } from './seq.js';
import { newId, colorFor, normRange, featureSpan, featureSequence } from './model.js';
import { BUILTIN_LIBRARY } from './feature-library.js';

export { BUILTIN_LIBRARY };

const clean = (s) => String(s ?? '').toUpperCase().replace(/[^A-Z]/g, '');

// Accept a library from JSON (a user's own, or an imported one) and drop what cannot be used.
export function normalizeLibrary(list) {
  const out = [];
  for (const e of Array.isArray(list) ? list : []) {
    const kind = e?.kind === 'protein' ? 'protein' : 'dna';
    const seq = clean(e?.seq);
    const minLen = kind === 'protein' ? 4 : 12;
    if (!e?.name || seq.length < minLen) continue;
    if (kind === 'dna' && /[^ACGTUN]/.test(seq)) continue;
    out.push({ name: String(e.name), type: e.type || 'misc_feature', kind, seq: kind === 'dna' ? seq.replace(/U/g, 'T') : seq, note: e.note ?? '' });
  }
  return out;
}

// Turn the annotated features of other documents into library entries, so a
// reader's own well-annotated plasmids teach the scanner about long elements.
export function libraryFromDocs(docs, { minLength = 18, maxLength = 6000 } = {}) {
  const seen = new Set();
  const out = [];
  for (const d of docs) {
    for (const f of d.features ?? []) {
      if (['source', 'gene', 'variation', 'crispr_site', 'primer_bind'].includes(f.type)) continue;
      let seq;
      try { seq = featureSequence(d, f).toUpperCase(); } catch { continue; }
      if (seq.length < minLength || seq.length > maxLength || /[^ACGT]/.test(seq) || seen.has(seq)) continue;
      seen.add(seq);
      out.push({ name: f.name, type: f.type, kind: 'dna', seq, note: `From ${d.name}`, fromDoc: d.id });
    }
  }
  return out;
}

// All places `query` matches `hay` within `budget` mismatches (hay positions of the match start).
function scanDna(hay, query, budget, limit) {
  const n = query.length;
  const out = [];
  if (budget === 0) {
    for (let at = hay.indexOf(query); at >= 0 && at < limit; at = hay.indexOf(query, at + 1)) out.push({ pos: at, mm: 0 });
    return out;
  }
  for (let pos = 0; pos < limit && pos + n <= hay.length; pos++) {
    let mm = 0;
    for (let k = 0; k < n; k++) {
      if (hay[pos + k] !== query[k] && ++mm > budget) break;
    }
    if (mm <= budget) out.push({ pos, mm });
  }
  return out;
}

function detectDna(doc, entry, rate, hits) {
  const L = doc.seq.length;
  const q = entry.seq;
  const n = q.length;
  if (n > L && !doc.circular) return;
  const hay = doc.circular ? doc.seq + doc.seq.slice(0, Math.max(0, n - 1)) : doc.seq;
  const limit = doc.circular ? L : L - n + 1;
  const budget = Math.floor(n * rate);
  const rc = revcomp(q);
  const strands = rc === q ? [[q, 1]] : [[q, 1], [rc, -1]];
  for (const [query, strand] of strands) {
    for (const { pos, mm } of scanDna(hay, query, budget, limit)) {
      const r = normRange(pos, n, L);
      hits.push({ name: entry.name, type: entry.type, kind: 'dna', strand, start: r.start, end: r.end, length: n, mismatches: mm, identity: (100 * (n - mm)) / n, note: entry.note, source: entry.fromDoc ? 'document' : 'library' });
    }
  }
}

function detectProtein(doc, entry, hits) {
  const L = doc.seq.length;
  const aaSeq = entry.seq;
  const need = aaSeq.length * 3;
  if (need > L && !doc.circular) return;
  for (const strand of [1, -1]) {
    const s = strand === 1 ? doc.seq : revcomp(doc.seq);
    const hay = doc.circular ? s + s.slice(0, need - 1) : s;
    for (let frame = 0; frame < 3; frame++) {
      const prot = translate(hay.slice(frame));
      for (let at = prot.indexOf(aaSeq); at >= 0; at = prot.indexOf(aaSeq, at + 1)) {
        const from = frame + at * 3; // in the (possibly reverse-complemented) reading strand
        if (from >= L) continue;
        const start = strand === 1 ? from : ((L - from - need) % L + L) % L;
        const r = normRange(start, need, L);
        hits.push({ name: entry.name, type: entry.type, kind: 'protein', strand, start: r.start, end: r.end, length: need, mismatches: 0, identity: 100, note: entry.note, source: 'library' });
      }
    }
  }
}

const sameSpan = (a, b) => a.start === b.start && a.end === b.end && a.strand === b.strand;

// Does hit `a` lie wholly inside hit `b` (same strand)? Positions are on the top strand.
function inside(a, b, L) {
  if (a.strand !== b.strand || a.length >= b.length) return false;
  const off = (a.start - b.start + L) % L;
  return off + a.length <= b.length;
}

export function detectFeatures(doc, library = BUILTIN_LIBRARY, { maxMismatchRate = 0.05, skipExisting = true } = {}) {
  const L = doc.seq.length;
  if (!L) return [];
  const raw = [];
  for (const entry of library) {
    if (entry.kind === 'protein') detectProtein(doc, entry, raw);
    else detectDna(doc, entry, maxMismatchRate, raw);
  }
  // the best-matching name for each span
  raw.sort((a, b) => b.identity - a.identity || b.length - a.length);
  const unique = [];
  for (const h of raw) if (!unique.some((u) => sameSpan(u, h))) unique.push(h);
  // a short element inside a longer one that was also found is part of it
  let hits = unique.filter((h) => !unique.some((o) => o !== h && inside(h, o, L)));
  if (skipExisting) {
    const have = doc.features.map((f) => ({ name: f.name.toLowerCase(), span: featureSpan(f), len: rangeLength(featureSpan(f), L) }));
    hits = hits.filter((h) => !have.some((f) => {
      if (f.name === h.name.toLowerCase()) return true; // already annotated, wherever
      const overlap = (Math.min(f.span.end, h.end) - Math.max(f.span.start, h.start));
      return f.span.start === h.start && f.span.end === h.end && overlap > 0;
    }));
  }
  return hits.sort((a, b) => a.start - b.start || b.length - a.length);
}

export function hitToFeature(hit) {
  const text = `Detected: ${hit.name}${hit.mismatches ? `, ${hit.mismatches} mismatch${hit.mismatches === 1 ? '' : 'es'}` : ''}${hit.note ? ` · ${hit.note}` : ''}`;
  return {
    id: newId(),
    name: hit.name,
    type: hit.type,
    strand: hit.strand,
    color: colorFor(hit.type),
    segments: [{ start: hit.start, end: hit.end }],
    qualifiers: [{ name: 'note', values: [{ kind: 'text', value: text }] }],
    extra: {},
  };
}
