// BublaGene document model.
//
// doc = {
//   id, name, seq (uppercase string), circular, flags (format flag byte),
//   features: [{ id, name, type, strand: 1|-1|0, color,
//                segments: [{start, end, color?, name?, extra?}], qualifiers: [{name, values:[{kind, value}]}], extra }],
//   primers:  [{ id, name, seq, description, extra }],
//   notes:    { description, organism, comments, created, modified, type, ... },
//   history:  [{ when, action, detail }],
//   source:   { format, raw?, ... }  (round-trip data for .dna)
//   seqEdited: bool
// }

import { mod, rangeLength, translate, proteinMass } from './seq.js';
import { isKnownTable, DEFAULT_TABLE } from './codon-tables.js';

let nextId = 1;
export const newId = () => `b${Date.now().toString(36)}${(nextId++).toString(36)}`;

// Feature keys BublaGene understands, with a display label and a small inline
// SVG icon path (24×24 viewBox). Unknown keys are still accepted and shown with
// a generic icon; nothing is discarded.
export const FEATURE_CATALOG = {
  CDS: { label: 'Coding sequence', icon: 'M4 6h16v12H4z M7 9h10M7 12h10M7 15h6', color: '#e0884c' },
  gene: { label: 'Gene', icon: 'M5 5h14M5 19h14M9 5v14M15 5v14', color: '#d6a14a' },
  promoter: { label: 'Promoter', icon: 'M4 5v14M7 5l8 7-8 7', color: '#3f9a7d' },
  terminator: { label: 'Terminator', icon: 'M5 6v12M9 6v12M13 8l7 4-7 4z', color: '#9b5fb8' },
  rep_origin: { label: 'Origin of replication', icon: 'M12 3a9 9 0 1 0 0 18a9 9 0 0 0 0-18z M9 12h6M12 9v6', color: '#d8c14a' },
  primer_bind: { label: 'Primer binding site', icon: 'M4 10h11M4 14h11M18 8v8', color: '#5b8fd6' },
  protein_bind: { label: 'Protein binding site', icon: 'M5 7h14v10H5z M8 11h8M8 14h5', color: '#4fb3b0' },
  RBS: { label: 'Ribosome binding site', icon: 'M5 8h14M5 16h14M8 8v8', color: '#7cb85b' },
  enhancer: { label: 'Enhancer', icon: 'M4 12h4l2-6 4 12 2-6h4', color: '#69a6c9' },
  'polyA_signal': { label: 'Poly(A) signal', icon: 'M5 8h14M8 12h11M11 16h8', color: '#b07aa1' },
  sig_peptide: { label: 'Signal peptide', icon: 'M4 16c4-8 12-8 16 0M7 16h4', color: '#e3a26f' },
  misc_feature: { label: 'Miscellaneous feature', icon: 'M5 5h14v14H5z', color: '#8c98a8' },
  LTR: { label: 'Long terminal repeat', icon: 'M7 5h10v14H7z M10 8h4M10 12h4M10 16h4', color: '#a6826a' },
  regulatory: { label: 'Regulatory region', icon: 'M4 12h3l2-5 3 10 3-10 2 5h3', color: '#5aa58f' },
  mRNA: { label: 'mRNA', icon: 'M5 7h14M5 12h14M5 17h9', color: '#c46f6f' },
  exon: { label: 'Exon', icon: 'M5 6h14v12H5z M11 6v12', color: '#c98b5a' },
  intron: { label: 'Intron', icon: 'M3 12h6M15 12h6M9 8c2 2 4 2 6 0M9 16c2-2 4-2 6 0', color: '#b4b4b4' },
  misc_RNA: { label: 'Miscellaneous RNA', icon: 'M5 7h14M5 12h14M5 17h9', color: '#c79a63' },
  ncRNA: { label: 'Non-coding RNA', icon: 'M12 4a8 8 0 1 0 0 16 M12 8v8M8 12h8', color: '#b98ac9' },
  tRNA: { label: 'Transfer RNA', icon: 'M6 5h12v6l-6 8-6-8z', color: '#7aa86a' },
  rRNA: { label: 'Ribosomal RNA', icon: 'M6 5h5v14H6z M13 5h5v14h-5z', color: '#6a9a8a' },
  mobile_element: { label: 'Mobile element', icon: 'M4 8h12v8H4z M16 10l4-2v8l-4-2', color: '#a98a6a' },
  source: { label: 'Source', icon: 'M5 5h14v14H5z M9 9h6v6H9z', color: '#8c98a8' },
  '-35_signal': { label: '−35 signal', icon: 'M5 8h14M5 16h10', color: '#5a9a8f' },
  '-10_signal': { label: '−10 signal', icon: 'M5 8h10M5 16h14', color: '#5a9a8f' },
  repeat_region: { label: 'Repeat region', icon: 'M5 9h14M5 12h14M5 15h14', color: '#9a9ad0' },
  'oriT': { label: 'Origin of transfer', icon: 'M4 12h16M15 8l5 4-5 4', color: '#d8c14a' },
  '5\'UTR': { label: "5' UTR", icon: 'M5 6h14M5 10h14M5 14h9', color: '#b0a0c0' },
  '3\'UTR': { label: "3' UTR", icon: 'M5 10h9M5 14h14M5 18h14', color: '#b0a0c0' },
};

export const FEATURE_TYPES = Object.keys(FEATURE_CATALOG);

const PALETTE = Object.fromEntries(Object.entries(FEATURE_CATALOG).map(([k, v]) => [k, v.color]));

export function colorFor(type) {
  return PALETTE[type] ?? '#8c98a8';
}

export function featureIcon(type) {
  return FEATURE_CATALOG[type]?.icon ?? FEATURE_CATALOG.misc_feature.icon;
}

export function featureLabel(type) {
  return FEATURE_CATALOG[type]?.label ?? String(type ?? 'feature');
}

// Default qualifiers for a newly created feature of this type.
export function defaultQualifiers(type) {
  if (type === 'CDS') return [{ name: 'transl_table', values: [{ kind: 'int', value: '1' }] }];
  return [];
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

// Is this segment part of the coding sequence? The format marks introns and
// other skipped stretches with type="gap"; everything else counts.
export const isCodingSegment = (s) => s.extra?.type !== 'gap';

// Sequence of a feature read 5'->3' in its own direction, over the whole span
// (gap segments included).
export function featureSequence(doc, f) {
  return joinSegments(doc, f, f.segments);
}

// How many bases to skip before the first codon. The format stores this as
// `readingFrame` on the feature (±1, ±2, ±3; the sign only marks the strand),
// matching the GenBank `codon_start` qualifier. Verified against the stored
// translations of 878 CDS features in the sample corpus.
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
  recomputeFeatures(doc);
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
  recomputeFeatures(doc);
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

export function qualifierInt(f, name, fallback = null) {
  const text = qualifierText(f, name);
  if (text === '') return fallback;
  const n = Number(text);
  return Number.isFinite(n) ? n : fallback;
}

// Genetic code used by a feature: its transl_table qualifier, else the
// document default, else NCBI table 1.
export function featureTableId(f, doc = null) {
  const q = qualifierInt(f, 'transl_table', null);
  if (q != null && isKnownTable(q)) return q;
  const e = Number(f?.extra?.transl_table);
  if (Number.isFinite(e) && isKnownTable(e)) return e;
  const d = Number(doc?.transl_table);
  if (Number.isFinite(d) && isKnownTable(d)) return d;
  return DEFAULT_TABLE;
}

// Reading-frame offset in bases (0, 1 or 2) from codon_start / stored frame.
export function phaseOffset(f) {
  return frameOffset(f) % 3;
}

// Translate a feature with its own genetic code. Partial codons are dropped.
export function translateFeature(doc, f, { initiator = false, stopSymbol = '*' } = {}) {
  const table = featureTableId(f, doc);
  return translate(codingSequence(doc, f), { table, initiator, stopSymbol });
}

// The translation a GenBank record normally stores: no initiator override and no
// terminal stop. Internal stops are kept (they are a fact about the sequence).
export function storedTranslation(doc, f) {
  const p = translateFeature(doc, f);
  return p.endsWith('*') ? p.slice(0, -1) : p;
}

// Recompute feature-derived qualifiers after an edit. Only touches values the
// feature already carries: a CDS gets /translation, a feature with /gc_content
// gets that too. Unknown qualifiers and their order are preserved.
export function recomputeFeature(doc, f) {
  if (!f) return f;
  if (f.type === 'CDS' || (f.type === 'gene' && (f.qualifiers ?? []).some((q) => q.name === 'translation'))) {
    const p = storedTranslation(doc, f);
    if ((f.qualifiers ?? []).some((q) => q.name === 'translation') || f.type === 'CDS') setQualifier(f, 'translation', p);
  }
  if ((f.qualifiers ?? []).some((q) => q.name === 'gc_content')) {
    const s = featureSequence(doc, f);
    setQualifier(f, 'gc_content', String(Math.round(gcContentLocal(s) * 10) / 10));
  }
  return f;
}

function gcContentLocal(s) {
  if (!s.length) return 0;
  let gc = 0;
  for (const ch of s.toUpperCase()) if (ch === 'G' || ch === 'C' || ch === 'S') gc++;
  return (100 * gc) / s.length;
}

// Recompute every derived qualifier in a document.
export function recomputeFeatures(doc) {
  for (const f of doc.features ?? []) recomputeFeature(doc, f);
  return doc;
}

// Protein molecular weight from a feature (or from a residue string).
export function proteinMolecularWeight(protein) {
  return proteinMass(protein);
}

export function featureProteinMass(doc, f) {
  return proteinMass(storedTranslation(doc, f));
}

export function featureGcContent(doc, f) {
  return gcContentLocal(featureSequence(doc, f));
}

// Is this segment's start/end uncertain (fuzzy location `<1` / `>50`)?
export function fuzzyStart(seg) {
  return Boolean(seg?.partialStart ?? seg?.extra?.partialStart);
}

export function fuzzyEnd(seg) {
  return Boolean(seg?.partialEnd ?? seg?.extra?.partialEnd);
}

// Merge segments that are adjacent on the top strand (no gap between them) while
// keeping every other property. Used after edits before saving.
export function normalizeSegments(segments) {
  const out = [];
  for (const s of segments ?? []) {
    const prev = out[out.length - 1];
    if (prev && prev.end === s.start && !fuzzyEnd(prev) && !fuzzyStart(s)) {
      out[out.length - 1] = { ...prev, end: s.end, partialEnd: s.partialEnd };
    } else {
      out.push({ ...s });
    }
  }
  return out;
}
