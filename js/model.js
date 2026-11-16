// BublaGene document model.
//
// doc = {
//   id, name, seq (uppercase string), circular, flags (SnapGene flag byte),
//   features: [{ id, name, type, strand: 1|-1|0, color,
//                segments: [{start, end, color?, name?, extra?}], qualifiers: [{name, values:[{kind, value}]}], extra }],
//   primers:  [{ id, name, seq, description, extra }],
//   notes:    { description, organism, comments, created, modified, type, ... },
//   history:  [{ when, action, detail }],
//   source:   { format, raw?, ... }  (round-trip data for .dna)
//   seqEdited: bool
// }

import { mod, rangeLength } from './seq.js';

let nextId = 1;
export const newId = () => `b${Date.now().toString(36)}${(nextId++).toString(36)}`;

export const FEATURE_TYPES = [
  'CDS', 'gene', 'promoter', 'terminator', 'rep_origin', 'primer_bind', 'misc_feature',
  'protein_bind', 'RBS', 'enhancer', 'polyA_signal', 'sig_peptide', 'mRNA', 'exon', 'intron',
  'LTR', 'regulatory', 'misc_RNA', 'ncRNA', 'tRNA', 'rRNA', 'mobile_element', 'source',
];

const PALETTE = {
  CDS: '#e0884c', gene: '#d6a14a', promoter: '#3f9a7d', terminator: '#9b5fb8',
  rep_origin: '#d8c14a', primer_bind: '#5b8fd6', protein_bind: '#4fb3b0', RBS: '#7cb85b',
  enhancer: '#69a6c9', polyA_signal: '#b07aa1', sig_peptide: '#e3a26f', misc_feature: '#8c98a8',
  LTR: '#a6826a', regulatory: '#5aa58f', mRNA: '#c46f6f', exon: '#c98b5a', intron: '#b4b4b4',
};

export function colorFor(type) {
  return PALETTE[type] ?? '#8c98a8';
}

export function newDocument({ name = 'Untitled', seq = '', circular = false } = {}) {
  const now = new Date();
  return {
    id: newId(),
    name,
    seq: seq.toUpperCase(),
    circular,
    flags: circular ? 0x1f : 0x1e,
    features: [],
    primers: [],
    notes: { type: 'Synthetic', description: '', organism: '', comments: '', created: now.toISOString(), modified: now.toISOString() },
    history: [{ when: now.toISOString(), action: 'created', detail: `New ${circular ? 'circular' : 'linear'} sequence, ${seq.length} bp` }],
    source: { format: 'new' },
    seqEdited: true,
  };
}

export function cloneDoc(doc) {
  const copy = structuredClone({ ...doc, source: undefined });
  copy.source = doc.source; // raw byte buffers are shared, never mutated
  return copy;
}

export function addHistory(doc, action, detail) {
  doc.history = doc.history || [];
  doc.history.push({ when: new Date().toISOString(), action, detail });
  doc.notes.modified = new Date().toISOString();
}

// Segments are always stored in top-strand order (as in .dna and GenBank files).
// Overall span of a feature: first segment start to last segment end.
export function featureSpan(f) {
  const segs = f.segments;
  if (!segs.length) return { start: 0, end: 0 };
  return { start: segs[0].start, end: segs[segs.length - 1].end };
}

export function featureLength(f, L) {
  return f.segments.reduce((n, s) => n + rangeLength(s, L), 0);
}

// Is this segment part of the coding sequence? SnapGene marks introns and other
// skipped stretches with type="gap"; everything else counts.
export const isCodingSegment = (s) => s.extra?.type !== 'gap';

// Sequence of a feature read 5'->3' in its own direction, over the whole span
// (gap segments included).
export function featureSequence(doc, f) {
  return joinSegments(doc, f, f.segments);
}

// How many bases to skip before the first codon. SnapGene stores this as
// `readingFrame` on the feature (±1, ±2, ±3; the sign only marks the strand),
// matching the GenBank `codon_start` qualifier. Verified against the stored
// translations of 878 CDS features in the SnapGene 9 sample files.
export function frameOffset(f) {
  const rf = Number(f.extra?.readingFrame);
  if (Number.isFinite(rf) && rf) return Math.abs(rf) - 1;
  const cs = Number(qualifierText(f, 'codon_start'));
  return Number.isFinite(cs) && cs ? cs - 1 : 0;
}

// The spliced coding sequence: gap segments (introns) left out and the reading
// frame applied. Use this for translation; featureSequence keeps everything.
export function codingSequence(doc, f) {
  return joinSegments(doc, f, f.segments.filter(isCodingSegment)).slice(frameOffset(f));
}

function joinSegments(doc, f, segs) {
  const top = segs
    .map((s) => (s.end > s.start ? doc.seq.slice(s.start, s.end) : doc.seq.slice(s.start) + doc.seq.slice(0, s.end)))
    .join('');
  return f.strand === -1 ? revcompLocal(top) : top;
}

function revcompLocal(s) {
  const C = { A: 'T', T: 'A', G: 'C', C: 'G', N: 'N', R: 'Y', Y: 'R', S: 'S', W: 'W', K: 'M', M: 'K', B: 'V', V: 'B', D: 'H', H: 'D' };
  let o = '';
  for (let i = s.length - 1; i >= 0; i--) o += C[s[i]] ?? s[i];
  return o;
}

// ---- editing ----------------------------------------------------------------

// Replace the bases in [start, end) (may wrap on circular) with `insert`.
// Features are shifted, trimmed or dropped so that they keep pointing at the same bases.
export function replaceRange(doc, start, end, insert) {
  const L = doc.seq.length;
  insert = insert.toUpperCase();
  let s = start;
  let e = end;
  // normalise a wrapping deletion by rotating the sequence first
  if (doc.circular && e < s) {
    rotateOrigin(doc, s);
    e = mod(e - s, L);
    s = 0;
  }
  const delLen = e - s;
  const shift = insert.length - delLen;
  doc.seq = doc.seq.slice(0, s) + insert + doc.seq.slice(e);
  const newL = doc.seq.length;
  // A feature start inside the replaced bases moves to just after the new text;
  // a feature end inside them is pulled back to the edit point.
  const adjStart = (p) => (p < s ? p : p >= e ? p + shift : s + insert.length);
  const adjEnd = (q) => (q <= s ? q : q > e ? q + shift : s);
  const keep = [];
  for (const f of doc.features) {
    const segs = [];
    for (const seg of f.segments) {
      const wraps = doc.circular && seg.end <= seg.start;
      const ns = adjStart(seg.start);
      const ne = adjEnd(seg.end);
      if (wraps) {
        const len = mod(ne - ns, newL) || newL;
        segs.push({ ...seg, ...normRange(ns, len, newL) });
      } else if (ne > ns) {
        segs.push({ ...seg, start: ns, end: ne });
      }
    }
    if (segs.length) keep.push({ ...f, segments: segs });
  }
  doc.features = keep;
  doc.seqEdited = true;
}

// Range of `len` bases starting at `start`. Ends are kept in 1..L, so a range
// ending at the origin has end === L, and a wrapping range has end <= start.
export function normRange(start, len, L) {
  start = mod(start, L);
  let end = start + len;
  if (end > L) end -= L;
  return { start, end };
}

// Make position p the new base 1 of a circular sequence.
export function rotateOrigin(doc, p) {
  const L = doc.seq.length;
  if (!doc.circular || !L) return;
  p = mod(p, L);
  if (!p) return;
  doc.seq = doc.seq.slice(p) + doc.seq.slice(0, p);
  for (const f of doc.features) {
    f.segments = f.segments.map((seg) => ({ ...seg, ...normRange(seg.start - p, rangeLength(seg, L), L) }));
  }
  doc.seqEdited = true;
}

export function reverseComplementDoc(doc) {
  const L = doc.seq.length;
  doc.seq = revcompLocal(doc.seq);
  for (const f of doc.features) {
    f.segments = f.segments
      .map((seg) => ({ ...seg, ...normRange(L - (seg.start + rangeLength(seg, L)), rangeLength(seg, L), L) }))
      .reverse();
    if (f.strand) f.strand = -f.strand;
  }
  doc.seqEdited = true;
}

// Features of `doc` lying entirely inside [from, from+len), shifted to start at
// `offset` in a new sequence (used for PCR products and region extraction).
export function featuresInRegion(doc, from, len, offset = 0) {
  const L = doc.seq.length;
  const out = [];
  for (const f of doc.features) {
    const segs = [];
    let ok = true;
    for (const s of f.segments) {
      const sl = rangeLength(s, L);
      const rel = doc.circular ? mod(s.start - from, L) : s.start - from;
      if (rel < 0 || rel + sl > len) { ok = false; break; }
      segs.push({ ...s, start: rel + offset, end: rel + offset + sl });
    }
    if (ok && segs.length) out.push({ ...structuredClone(f), id: newId(), segments: segs });
  }
  return out;
}

export function qualifierText(f, name) {
  const q = f.qualifiers?.find((x) => x.name === name);
  if (!q) return '';
  return q.values.map((v) => v.value).join('; ');
}

export function setQualifier(f, name, text) {
  f.qualifiers = (f.qualifiers || []).filter((q) => q.name !== name);
  if (text) f.qualifiers.push({ name, values: [{ kind: 'text', value: text }] });
}
