// Search and filtering across sequence, features and qualifiers.

import { findPattern, iupacRegexSource, mod } from './seq.js';

function qualText(f) {
  return (f.qualifiers ?? []).map((q) => `${q.name}=${q.values.map((v) => v.value).join(';')}`).join(' ');
}

export function searchSequence(doc, query, {
  bothStrands = true, circular = doc.circular, regex = false, caseSensitive = false, max = 500,
} = {}) {
  if (!query) return [];
  const seq = caseSensitive ? doc.seq : doc.seq.toUpperCase();
  const needle = caseSensitive ? query : query.toUpperCase();
  if (regex) {
    const re = new RegExp(needle, caseSensitive ? 'g' : 'gi');
    const out = [];
    let m;
    const hay = circular ? seq + seq.slice(0, Math.max(0, needle.length)) : seq;
    while ((m = re.exec(hay)) && out.length < max) {
      if (m.index >= seq.length) break;
      out.push({ kind: 'sequence', strand: 1, start: m.index, end: circular ? mod(m.index + m[0].length, seq.length) || m.index + m[0].length : m.index + m[0].length, text: m[0] });
      if (m[0].length === 0) re.lastIndex++;
    }
    return out;
  }
  let hits;
  try {
    hits = findPattern(seq, needle, circular);
  } catch {
    return []; // not a nucleotide pattern; feature/qualifier search still runs
  }
  hits = hits.filter((h) => bothStrands || h.strand === 1);
  return hits.slice(0, max).map((h) => ({ kind: 'sequence', strand: h.strand, start: h.start, end: h.end, text: needle }));
}

export function searchFeatures(doc, query, { caseSensitive = false, max = 500 } = {}) {
  if (!query) return [];
  const q = caseSensitive ? query : query.toLowerCase();
  const matches = (text) => (caseSensitive ? text : text.toLowerCase()).includes(q);
  const out = [];
  for (const f of doc.features ?? []) {
    if (matches(f.name)) out.push({ kind: 'feature', field: 'name', featureId: f.id, name: f.name, type: f.type, text: f.name });
    if (matches(f.type)) out.push({ kind: 'feature', field: 'type', featureId: f.id, name: f.name, type: f.type, text: f.type });
    for (const qual of f.qualifiers ?? []) {
      for (const value of qual.values ?? []) {
        if (matches(String(value.value ?? '')) || matches(qual.name)) {
          out.push({ kind: 'qualifier', field: qual.name, featureId: f.id, name: f.name, type: f.type, text: String(value.value ?? '') });
        }
      }
    }
    if (out.length >= max) break;
  }
  return out.slice(0, max);
}

export function filterFeatures(doc, { type = null, strand = null, text = null, hasQualifier = null, region = null } = {}) {
  return (doc.features ?? []).filter((f) => {
    if (type && f.type !== type) return false;
    if (strand != null && f.strand !== Number(strand)) return false;
    if (hasQualifier && !(f.qualifiers ?? []).some((q) => q.name === hasQualifier)) return false;
    if (text && !`${f.name} ${f.type} ${qualText(f)}`.toLowerCase().includes(String(text).toLowerCase())) return false;
    if (region) {
      const spans = f.segments ?? [];
      if (!spans.some((s) => s.end > region.start && s.start < region.end)) return false;
    }
    return true;
  });
}

export function searchAll(doc, query, options = {}) {
  return [...searchSequence(doc, query, options), ...searchFeatures(doc, query, options)];
}

export function featureTypes(doc) {
  return [...new Set((doc.features ?? []).map((f) => f.type))].sort();
}

export function qualifierKeys(doc) {
  return [...new Set((doc.features ?? []).flatMap((f) => (f.qualifiers ?? []).map((q) => q.name)))].sort();
}

export { iupacRegexSource };
