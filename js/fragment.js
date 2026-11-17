// Fragment copy/paste and feature reordering.

import { featuresInRegion, newId, featureSpan, normRange, replaceRange } from './model.js';
import { rangeLength, subseq, mod } from './seq.js';

// Capture the bases and the features fully inside a selection.
export function extractFragment(doc, range, name = null) {
  if (!range) return null;
  const len = rangeLength(range, doc.seq.length);
  if (!len) return null;
  const seq = subseq(doc.seq, range.start, range.end, doc.circular).toUpperCase();
  const features = featuresInRegion(doc, range.start, len, 0).map((f) => ({ ...f, id: newId() }));
  return {
    name: name ?? `${doc.name} ${range.start + 1}–${(range.start + len - 1) % Math.max(1, doc.seq.length) + 1}`,
    seq,
    features,
    circular: false,
    length: seq.length,
  };
}

// Insert a captured fragment at `at`, replacing nothing. Existing features shift
// with the sequence edit; pasted features are added at the insertion point.
export function pasteFragment(doc, at, fragment) {
  if (!fragment?.seq) throw new Error('Nothing to paste');
  const L = doc.seq.length;
  const position = doc.circular ? mod(at, L) : Math.max(0, Math.min(L, at));
  const before = doc.seq.length;
  replaceRange(doc, position, position, fragment.seq);
  const newL = doc.seq.length;
  for (const f of fragment.features ?? []) {
    const segments = (f.segments ?? []).map((s) => {
      const start = position + s.start;
      const end = position + s.end;
      return doc.circular ? { ...s, ...normRange(start, end - s.start, newL) } : { ...s, start, end };
    }).filter((s) => s.end > s.start);
    if (segments.length) doc.features.push({ ...structuredClone(f), id: newId(), segments });
  }
  return { position, inserted: fragment.seq, length: doc.seq.length - before };
}

// Move a feature so it sits before `beforeId` in the feature list.
export function reorderFeature(doc, id, beforeId = null) {
  const from = doc.features.findIndex((f) => f.id === id);
  if (from < 0) return false;
  const [feature] = doc.features.splice(from, 1);
  let to = beforeId ? doc.features.findIndex((f) => f.id === beforeId) : -1;
  if (to < 0) to = doc.features.length;
  doc.features.splice(to, 0, feature);
  return true;
}

// Selected feature spans, merged when they touch, for copying or deleting.
export function selectedSpans(doc, features) {
  const spans = features.map((f) => featureSpan(f)).filter((s) => s.end > s.start).sort((a, b) => a.start - b.start);
  const out = [];
  for (const s of spans) {
    const last = out[out.length - 1];
    if (last && s.start <= last.end) last.end = Math.max(last.end, s.end);
    else out.push({ ...s });
  }
  return out;
}
