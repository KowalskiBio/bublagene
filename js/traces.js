// Trace trimming and conversion, shared by ABI, SCF and ZTR readers.
//
// Mott trimming uses a sliding-window mean Phred score; Phred trimming is the
// simpler first/last base above a threshold. Both return base coordinates and
// can add trim features to a document.

import { readAB1 } from './formats/ab1.js';
import { readSCF } from './formats/scf.js';
import { readZtr } from './formats/ztr.js';
import { newDocument, newId } from './model.js';

const ZTR_MAGIC = [0xae, 0x5a, 0x54, 0x52, 0x0d, 0x0a, 0x1a, 0x0a];

export function traceFormat(buffer) {
  const b = buffer instanceof Uint8Array ? buffer : new Uint8Array(buffer);
  if (String.fromCharCode(b[0], b[1], b[2], b[3]) === 'ABIF') return 'AB1';
  if (String.fromCharCode(b[0], b[1], b[2], b[3]) === '.scf') return 'SCF';
  if (ZTR_MAGIC.every((x, i) => b[i] === x)) return 'ZTR';
  return null;
}

export async function readTrace(buffer, fileName = 'trace') {
  const format = traceFormat(buffer);
  if (format === 'AB1') return readAB1(buffer, fileName);
  if (format === 'SCF') return readSCF(buffer, fileName);
  if (format === 'ZTR') return readZtr(buffer instanceof Uint8Array ? buffer : new Uint8Array(buffer));
  throw new Error('Unrecognised trace format (expected ABI, SCF or ZTR)');
}

export function qualitiesOf(trace) {
  const q = trace.confidence ?? trace.qualities;
  if (q && q.length >= trace.bases.length) return Array.from(q).slice(0, trace.bases.length);
  return Array.from({ length: trace.bases.length }, () => 20);
}

export function qualityStats(trace) {
  const q = qualitiesOf(trace);
  if (!q.length) return { count: 0, min: 0, max: 0, mean: 0, q20: 0, q30: 0 };
  const sorted = [...q].sort((a, b) => a - b);
  return {
    count: q.length,
    min: sorted[0],
    max: sorted[sorted.length - 1],
    mean: q.reduce((n, x) => n + x, 0) / q.length,
    q20: sorted[Math.floor(sorted.length * 0.2)],
    q30: sorted[Math.floor(sorted.length * 0.3)],
  };
}

// Return {start, end} base coordinates. With `mott`, windows whose mean Phred
// is below `cutoff` are removed from each end.
export function trimTrace(trace, { mode = 'mott', window = 10, cutoff = 20 } = {}) {
  const q = qualitiesOf(trace);
  const n = q.length;
  if (!n) return { start: 0, end: 0, length: 0, removed: { start: 0, end: 0 }, mode, cutoff, window, mean: 0 };
  let start = 0;
  let end = n;
  if (mode === 'phred') {
    start = q.findIndex((x) => x >= cutoff);
    if (start < 0) start = n;
    end = n - [...q].reverse().findIndex((x) => x >= cutoff);
    if (end > n) end = n;
  } else {
    const w = Math.max(1, Math.min(window, n));
    const means = [];
    for (let i = 0; i + w <= n; i++) {
      let sum = 0;
      for (let k = 0; k < w; k++) sum += q[i + k];
      means.push(sum / w);
    }
    const first = means.findIndex((m) => m >= cutoff);
    start = first < 0 ? n : first;
    let last = -1;
    for (let i = means.length - 1; i >= 0; i--) if (means[i] >= cutoff) { last = i; break; }
    end = last < 0 ? 0 : last + w;
  }
  if (end < start) end = start;
  return {
    start,
    end,
    length: Math.max(0, end - start),
    removed: { start, end: n - end },
    mode,
    cutoff,
    window,
    mean: q.slice(start, end).reduce((a, b) => a + b, 0) / Math.max(1, end - start),
  };
}

// Add low-quality trim features to a document that carries the trace bases.
export function trimFeatures(doc, trim) {
  const L = doc.seq.length;
  const out = [];
  const add = (start, end, side) => {
    if (end <= start) return;
    out.push({
      id: newId(), name: `trimmed ${side}`, type: 'misc_feature', strand: 1, color: '#8c98a8',
      segments: [{ start, end, color: '#8c98a8' }],
      qualifiers: [{ name: 'note', values: [{ kind: 'text', value: `low quality ${side} (${trim.mode}, cutoff ${trim.cutoff})` }] }],
      extra: {},
    });
  };
  add(0, Math.min(trim.start, L), 'start');
  add(Math.min(trim.end, L), L, 'end');
  return out;
}

export function traceToDocument(trace, { name = null, trimmed = false } = {}) {
  const doc = newDocument({ name: name ?? trace.name ?? 'trace', seq: trace.bases.toUpperCase(), circular: false });
  doc.source = {
    format: 'trace', traceFormat: trace.format, fileName: trace.name,
    positions: [...(trace.positions ?? [])], qualities: qualitiesOf(trace),
  };
  doc.trace = {
    format: trace.format,
    name: trace.name,
    sampleCount: trace.sampleCount,
    qualities: qualitiesOf(trace),
    positions: [...(trace.positions ?? [])],
  };
  doc.notes.description = `${trace.format} trace${trace.name ? ` ${trace.name}` : ''}: ${trace.bases.length} bases, ${trace.sampleCount} samples`;
  if (trimmed) doc.features = trimFeatures(doc, trimTrace(trace, { mode: 'phred', cutoff: 20 }));
  return doc;
}

// Slice a trace to a base range (keeps sample coordinates and called bases).
export function sliceTrace(trace, start, end) {
  const positions = (trace.positions ?? []).slice(start, end);
  const sampleStart = positions.length ? positions[0] : 0;
  const sampleEnd = positions.length ? (positions[positions.length - 1] + 1) : 0;
  const samples = {};
  for (const base of ['A', 'C', 'G', 'T']) {
    const src = trace.samples?.[base];
    samples[base] = src ? src.slice(Math.max(0, sampleStart), Math.max(sampleStart, sampleEnd)) : new Uint16Array(0);
  }
  return {
    ...trace,
    bases: trace.bases.slice(start, end),
    positions: positions.map((p) => p - sampleStart),
    confidence: (trace.confidence ?? []).slice(start, end),
    qualities: (trace.qualities ?? []).slice(start, end),
    samples,
    sampleCount: Math.max(0, sampleEnd - sampleStart),
    sourceRange: { start, end },
  };
}

// Map a sample coordinate to the nearest called base, for click navigation.
export function baseAtSample(trace, sample) {
  const positions = trace.positions ?? [];
  if (!positions.length) return null;
  let best = 0;
  let distance = Infinity;
  for (let i = 0; i < positions.length; i++) {
    const d = Math.abs(positions[i] - sample);
    if (d < distance) { distance = d; best = i; }
  }
  return { index: best, base: trace.bases[best], position: positions[best], distance };
}
