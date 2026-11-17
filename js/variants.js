// Variant calling between a reference and a query sequence.
//
// Variants use reference-anchored, 0-based coordinates. Indels are normalised
// (common prefix/suffix trimmed, homopolymer indels shifted left) so the same
// event has one representation. `variantsToFeatures` / `featuresToVariants`
// round-trip the calls through ordinary document features.

import { alignGlobal } from './align.js';
import { mod } from './seq.js';
import { newId } from './model.js';

const TYPE_OF = (ref, alt) => (ref.length === 0 ? 'INS' : alt.length === 0 ? 'DEL' : ref.length === 1 && alt.length === 1 ? 'SNP' : 'MNV');

// Left-align and trim an indel. `reference` is the full reference sequence.
export function normalizeVariant(variant, reference) {
  let { position, ref, alt } = variant;
  if (!ref && !alt) return { ...variant, type: 'NOP' };
  // trim common suffix
  while (ref && alt && ref.at(-1) === alt.at(-1)) { ref = ref.slice(0, -1); alt = alt.slice(0, -1); }
  // trim common prefix
  while (ref && alt && ref[0] === alt[0]) { ref = ref.slice(1); alt = alt.slice(1); position++; }
  if (!ref && !alt) return { ...variant, position, ref, alt, type: 'NOP' };
  if (!ref || !alt) {
    let allele = ref || alt;
    let guard = 0;
    while (position > 0 && allele.length && reference[position - 1] === allele.at(-1) && guard++ < 10000) {
      position--;
      allele = allele.at(-1) + allele.slice(0, -1);
    }
    if (ref) ref = allele; else alt = allele;
  }
  return { ...variant, position, ref, alt, type: TYPE_OF(ref, alt) };
}

function parseAlignment(alignedA, alignedB, startA) {
  const out = [];
  let i = 0;
  let refPos = startA;
  while (i < alignedA.length) {
    const a = alignedA[i];
    const b = alignedB[i];
    if (a !== '-' && b !== '-' && a === b) { i++; refPos++; continue; }
    if (a === '-' && b === '-') { i++; continue; }
    if (a === '-') {
      // insertion relative to the reference
      let alt = '';
      const position = refPos;
      while (i < alignedA.length && alignedA[i] === '-' && alignedB[i] !== '-') { alt += alignedB[i]; i++; }
      out.push({ type: 'INS', position, ref: '', alt });
      continue;
    }
    if (b === '-') {
      let ref = '';
      const position = refPos;
      while (i < alignedA.length && alignedB[i] === '-' && alignedA[i] !== '-') { ref += alignedA[i]; i++; refPos++; }
      out.push({ type: 'DEL', position, ref, alt: '' });
      continue;
    }
    let ref = '';
    let alt = '';
    const position = refPos;
    while (i < alignedA.length && alignedA[i] !== '-' && alignedB[i] !== '-' && alignedA[i] !== alignedB[i]) {
      ref += alignedA[i];
      alt += alignedB[i];
      i++;
      refPos++;
    }
    out.push({ type: TYPE_OF(ref, alt), position, ref, alt });
  }
  return out;
}

export function callVariants(reference, query, { circular = false, band = null, ...alignOpts } = {}) {
  const ref = String(reference).toUpperCase();
  const qry = String(query).toUpperCase();
  const L = ref.length;
  if (!circular) {
    const alignment = alignGlobal(ref, qry, { band, ...alignOpts });
    const variants = parseAlignment(alignment.alignedA, alignment.alignedB, alignment.startA)
      .map((v) => normalizeVariant(v, ref))
      .filter((v) => v.type !== 'NOP');
    return { variants, alignment, score: alignment.score, identity: alignment.identity, circular };
  }
  // Try two rotations so a variant at the origin is still seen.
  const candidates = [];
  for (const shift of [0, Math.floor(L / 2)]) {
    const rotated = ref.slice(shift) + ref.slice(0, shift);
    const alignment = alignGlobal(rotated, qry, { band, ...alignOpts });
    const variants = parseAlignment(alignment.alignedA, alignment.alignedB, alignment.startA)
      .map((v) => normalizeVariant(v, rotated))
      .filter((v) => v.type !== 'NOP')
      .map((v) => ({ ...v, position: mod(v.position + shift, L) }));
    candidates.push({ variants, alignment, score: alignment.score, shift });
  }
  candidates.sort((a, b) => b.score - a.score);
  const best = candidates[0];
  return { variants: best.variants, alignment: best.alignment, score: best.score, identity: best.alignment.identity, circular, shift: best.shift };
}

export function variantKey(v) {
  return `${v.type}:${v.position}:${v.ref}:${v.alt}`;
}

function qualifier(name, value) {
  return { name, values: [{ kind: 'text', value: String(value) }] };
}

export function variantsToFeatures(doc, variants) {
  const L = doc.seq.length;
  return variants.map((v) => {
    const start = L ? Math.min(v.position, L - 1) : 0;
    const end = v.type === 'DEL'
      ? Math.min(L, start + v.ref.length)
      : v.type === 'INS'
        ? Math.min(L, start + 1)
        : Math.min(L, start + Math.max(1, v.ref.length));
    return {
      id: newId(),
      name: `${v.type} ${v.position + 1}`,
      type: 'variation',
      strand: 1,
      color: '#c06a6a',
      segments: [{ start, end: Math.max(end, start + 1), color: '#c06a6a' }],
      qualifiers: [
        qualifier('variant_type', v.type),
        qualifier('reference_base', v.ref),
        qualifier('variant_base', v.alt),
        qualifier('reference_position', v.position + 1),
      ],
      extra: {},
    };
  });
}

export function featuresToVariants(doc) {
  const out = [];
  for (const f of doc.features ?? []) {
    if (f.type !== 'variation') continue;
    const get = (name) => f.qualifiers?.find((q) => q.name === name)?.values?.[0]?.value ?? '';
    const type = get('variant_type');
    const ref = get('reference_base');
    const alt = get('variant_base');
    const pos = Number(get('reference_position'));
    if (!type) continue;
    out.push({ type, position: Number.isFinite(pos) && pos > 0 ? pos - 1 : f.segments[0]?.start ?? 0, ref, alt });
  }
  return out;
}

export function formatVariant(v) {
  const where = `${v.position + 1}`;
  if (v.type === 'SNP') return `${where} ${v.ref}>${v.alt}`;
  if (v.type === 'INS') return `${where} ins ${v.alt}`;
  if (v.type === 'DEL') return `${where} del ${v.ref}`;
  return `${where} ${v.ref}>${v.alt}`;
}

// Apply variants to a reference (for tests and for building a query).
export function applyVariants(reference, variants) {
  const ref = String(reference).toUpperCase();
  const sorted = [...variants].sort((a, b) => b.position - a.position || b.ref.length - a.ref.length);
  let out = ref;
  for (const v of sorted) {
    out = out.slice(0, v.position) + (v.alt ?? '') + out.slice(v.position + (v.ref?.length ?? 0));
  }
  return out;
}
