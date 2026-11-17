// Site-specific recombination (Gateway-style att sites, Cre/lox, Flp/FRT) and
// TOPO-style insertion.
//
// Recombination is modelled as the double crossover it is: two circles (or a
// circle and a line) that each carry two matching sites swap the stretch
// between them. The crossover happens inside the sites' shared core, so the new
// junction is [acceptor side of the core] + core + [insert side of the core].
// Nothing here depends on knowing what the arms of a site look like: the sites
// are taken from the documents' own annotations, and the core is whatever the
// two sites have in common.

import { revcomp, mod, rangeLength } from './seq.js';
import { newDocument, newId, addHistory, cloneDoc, reverseComplementDoc, rotateOrigin, featureSpan, replaceRange } from './model.js';
import { extractFragment, pasteFragment } from './fragment.js';
import { detectFeatures, BUILTIN_LIBRARY } from './detect.js';
import { CLONING_SCHEMA } from './cloning.js';

const SITE_NAME = /^(att[BPLR][1-5]|lox[\w]*|FRT[\w]*)$/i;
const ATT = /^att([BPLR])([1-5])$/i;

// The key two sites must share to recombine: the number of an att site, or the
// full name of a lox or FRT site.
export function siteKey(name) {
  const m = ATT.exec(name);
  return m ? `att${m[2]}` : String(name).toUpperCase();
}

// Sites in a document: annotated features first, then built-in att sites found by sequence.
export function findRecombinationSites(doc) {
  const L = doc.seq.length;
  const out = [];
  for (const f of doc.features) {
    if (!SITE_NAME.test(f.name) || !f.segments.length) continue;
    const span = featureSpan(f);
    const len = rangeLength(span, L);
    out.push({ name: f.name, key: siteKey(f.name), start: span.start, strand: f.strand === -1 ? -1 : 1, length: len, seq: slice(doc, span.start, len) });
  }
  const known = new Set(out.map((s) => s.name.toLowerCase()));
  const scan = BUILTIN_LIBRARY.filter((e) => ATT.test(e.name) && !known.has(e.name.toLowerCase()));
  for (const h of detectFeatures(doc, scan, { maxMismatchRate: 0, skipExisting: false })) {
    out.push({ name: h.name, key: siteKey(h.name), start: h.start, strand: h.strand, length: h.length, seq: slice(doc, h.start, h.length), detected: true });
  }
  return out.sort((a, b) => a.start - b.start);
}

function slice(doc, start, len) {
  const L = doc.seq.length;
  if (!L) return '';
  const s = mod(start, L);
  return len > L ? '' : (doc.seq + doc.seq).slice(s, s + len);
}

// Longest stretch two sequences share, as [index in a, index in b, length].
function longestCommon(a, b) {
  let best = [0, 0, 0];
  let prev = new Array(b.length + 1).fill(0);
  for (let i = 1; i <= a.length; i++) {
    const cur = new Array(b.length + 1).fill(0);
    for (let j = 1; j <= b.length; j++) {
      if (a[i - 1] === b[j - 1]) {
        cur[j] = prev[j - 1] + 1;
        if (cur[j] > best[2]) best = [i - cur[j], j - cur[j], cur[j]];
      }
    }
    prev = cur;
  }
  return best;
}

// Features of the stretch [start, start+len) of a (possibly circular) document,
// positioned relative to that stretch. Features cut by its edges are trimmed.
function featuresOfArc(doc, start, len) {
  const c = cloneDoc(doc);
  const L = c.seq.length;
  if (c.circular) rotateOrigin(c, start);
  const base = c.circular ? 0 : start;
  const out = [];
  for (const f of c.features) {
    const segs = [];
    for (const seg of f.segments) {
      const pieces = seg.end > seg.start ? [[seg.start, seg.end]] : [[seg.start, L], [0, seg.end]];
      for (const [a, b] of pieces) {
        const s0 = Math.max(a, base);
        const s1 = Math.min(b, base + len);
        if (s1 > s0) segs.push({ ...structuredClone(seg), start: s0 - base, end: s1 - base });
      }
    }
    if (segs.length) out.push({ ...structuredClone(f), id: newId(), segments: segs });
  }
  return out;
}

const arc = (doc, from, to) => {
  const L = doc.seq.length;
  const len = mod(to - from, L);
  return { start: mod(from, L), len, seq: len ? slice(doc, from, len) : '' };
};

// Orient a document so its two sites both read forward; reject a mixed pair.
function orient(doc, keys) {
  let work = doc;
  for (let pass = 0; pass < 2; pass++) {
    const sites = keys.map((k) => findRecombinationSites(work).find((s) => s.key === k));
    if (sites.some((s) => !s)) throw new Error(`${doc.name}: missing a ${keys.find((k, i) => !sites[i])} site`);
    if (sites[0].strand !== sites[1].strand) throw new Error(`${doc.name}: the two sites face opposite ways, which would invert the insert`);
    if (sites[0].strand === 1) return { doc: work, sites };
    work = cloneDoc(doc);
    reverseComplementDoc(work);
  }
  throw new Error(`${doc.name}: could not orient the sites`);
}

const NEW_SITE = { LR: 'attB', BP: 'attL', PB: 'attR', RL: 'attP' };
const kindOf = (name) => ATT.exec(name)?.[1]?.toUpperCase() ?? null;

/**
 * Swap the stretch between two matching sites of `insertDoc` into `acceptorDoc`.
 * Gateway LR: insertDoc is the entry clone, acceptorDoc the destination vector.
 * Gateway BP: insertDoc is the attB-flanked fragment, acceptorDoc the donor.
 */
export function recombine(insertDoc, acceptorDoc, { name = null, minCore = 10 } = {}) {
  if (!acceptorDoc.circular) throw new Error(`${acceptorDoc.name} must be circular to receive the insert`);
  const sa = findRecombinationSites(insertDoc);
  const sb = findRecombinationSites(acceptorDoc);
  const shared = [...new Set(sa.map((s) => s.key))].filter((k) => sb.some((s) => s.key === k)).sort();
  if (shared.length < 2) {
    throw new Error(shared.length
      ? `Only one site is shared (${shared[0]}); recombination needs two to move a segment`
      : 'No matching recombination sites: annotate att, lox or FRT sites on both files');
  }
  const keys = shared.slice(0, 2);
  const A = orient(insertDoc, keys);
  const B = orient(acceptorDoc, keys);

  // The crossover core of each pair: what the two copies of a site have in common.
  const cores = keys.map((k, i) => {
    const [ia, ib, n] = longestCommon(A.sites[i].seq, B.sites[i].seq);
    if (n < minCore) throw new Error(`${k}: the two sites share only ${n} bp, too little to recombine`);
    return { key: k, seq: A.sites[i].seq.slice(ia, ia + n), a: mod(A.sites[i].start + ia, A.doc.seq.length), b: mod(B.sites[i].start + ib, B.doc.seq.length), n };
  });
  const [c1, c2] = cores;

  const region = arc(A.doc, c1.a + c1.n, c2.a);          // insert: between the cores
  const backbone = arc(B.doc, c2.b + c2.n, c1.b);        // acceptor: the part kept
  if (!A.doc.circular && c2.a < c1.a) throw new Error(`${A.doc.name}: site ${keys[1]} comes before ${keys[0]}; flip the file first`);
  if (!region.len && !A.doc.circular) throw new Error('The two sites leave nothing between them');

  const seq = backbone.seq + c1.seq + region.seq + c2.seq;
  const kindA = kindOf(A.sites[0].name);
  const kindB = kindOf(B.sites[0].name);
  const pair = kindA && kindB ? `${kindA}${kindB}` : null;
  const productName = name ?? `${B.doc.name} + ${A.doc.name}`;
  const doc = newDocument({ name: productName, seq, circular: true });

  // features: the acceptor's backbone, then the insert's stretch, offset into place
  const place = (features, offset) => {
    for (const f of features) {
      f.segments = f.segments.map((s) => ({ ...s, start: s.start + offset, end: s.end + offset }));
      doc.features.push(f);
    }
  };
  place(featuresOfArc(B.doc, backbone.start, backbone.len).filter((f) => !SITE_NAME.test(f.name)), 0);
  place(featuresOfArc(A.doc, mod(c1.a + c1.n, A.doc.seq.length), region.len).filter((f) => !SITE_NAME.test(f.name)), backbone.len + c1.n);
  const siteName = (i) => {
    const letter = pair ? NEW_SITE[pair] : null;
    return letter ? `${letter}${keys[i].replace(/^att/i, '')}` : keys[i];
  };
  const junctionAt = [backbone.len, backbone.len + c1.n + region.len];
  cores.forEach((c, i) => {
    doc.features.push({
      id: newId(), name: siteName(i), type: 'misc_feature', strand: 1, color: '#8c98a8',
      segments: [{ start: junctionAt[i], end: junctionAt[i] + c.n }],
      qualifiers: [{ name: 'note', values: [{ kind: 'text', value: `Recombined core (${c.n} bp) of ${A.sites[i].name} × ${B.sites[i].name}` }] }], extra: {},
    });
  });

  // the other product, when the insert came from a circle
  let byproduct = null;
  if (A.doc.circular) {
    const aBack = arc(A.doc, c2.a + c2.n, c1.a);
    const bRegion = arc(B.doc, c1.b + c1.n, c2.b);
    byproduct = newDocument({ name: `${A.doc.name} + ${B.doc.name} (by-product)`, seq: aBack.seq + c1.seq + bRegion.seq + c2.seq, circular: true });
  }

  const created = new Date().toISOString();
  doc.cloning = {
    schema: CLONING_SCHEMA, type: 'recombination', operation: 'site-specific-recombination', created, circular: true,
    product: { name: doc.name, seqLen: doc.seq.length },
    inputs: [
      { name: B.doc.name, seqLen: B.doc.seq.length, role: 'backbone', orientation: B.doc === acceptorDoc ? 'forward' : 'reverse' },
      { name: A.doc.name, seqLen: A.doc.seq.length, role: 'insert', orientation: A.doc === insertDoc ? 'forward' : 'reverse' },
    ],
    junctions: cores.map((c, i) => ({
      position: junctionAt[i], type: 'recombination', overhang: c.seq, enzymes: [], leftInput: i === 0 ? 0 : 1, rightInput: i === 0 ? 1 : 0,
      sequence: { top: c.seq, bottom: revcomp(c.seq), start: junctionAt[i] }, regeneratedSites: [],
    })),
  };
  addHistory(doc, 'Recombination', `${A.doc.name} × ${B.doc.name} at ${keys.join(' and ')} → ${doc.seq.length} bp`);
  doc.notes.description = `Segment between ${keys.join(' and ')} of ${A.doc.name} moved into ${B.doc.name}.`;
  const warnings = [];
  if (A.doc !== insertDoc || B.doc !== acceptorDoc) warnings.push('A file was reverse-complemented so its sites read the same way.');
  if (shared.length > 2) warnings.push(`${shared.length} site pairs match; used ${keys.join(' and ')}.`);
  return { doc, byproduct, sites: { keys, cores }, warnings };
}

// ---- TOPO and other insertion cloning ---------------------------------------

const KINDS = ['TA', 'blunt', 'directional'];

/**
 * Put an insert into a vector at a point or over a range, the way TOPO-type
 * cloning does. TA and blunt cloning are not directional, so the orientation is
 * the caller's choice; directional cloning needs the insert to begin with CACC.
 */
export function topoClone(vector, insert, { kind = 'TA', range = null, orientation = 'forward', name = null } = {}) {
  if (!KINDS.includes(kind)) throw new Error(`Unknown cloning kind ${kind}`);
  if (!vector?.seq) throw new Error('Choose a vector');
  const L = vector.seq.length;
  const at = range ?? { start: 0, end: 0 };
  if (!range) throw new Error('Select where the insert goes in the vector');
  if (at.end < at.start) throw new Error('The cloning site cannot wrap around the origin');
  const ins = cloneDoc(insert);
  if (orientation === 'reverse') reverseComplementDoc(ins);
  if (!ins.seq) throw new Error('The insert is empty');
  if (kind === 'directional' && !ins.seq.startsWith('CACC')) {
    throw new Error(orientation === 'forward' && revcomp(insert.seq).startsWith('CACC')
      ? 'The CACC end is on the other strand: choose the reverse orientation or flip the insert'
      : 'Directional cloning needs the insert to begin with CACC');
  }
  const fragment = extractFragment(ins, { start: 0, end: ins.seq.length });
  const doc = cloneDoc(vector);
  replaceRange(doc, at.start, at.end, '');
  pasteFragment(doc, at.start, fragment);
  doc.id = newId();
  doc.name = name ?? `${vector.name} + ${insert.name}`;
  doc.source = { format: 'new' };
  doc.cloning = {
    schema: CLONING_SCHEMA, type: 'topo', operation: `${kind.toLowerCase()}-topo-cloning`, created: new Date().toISOString(), circular: !!vector.circular,
    product: { name: doc.name, seqLen: doc.seq.length },
    inputs: [
      { name: vector.name, seqLen: L, role: 'backbone', orientation: 'forward' },
      { name: insert.name, seqLen: insert.seq.length, role: 'insert', orientation },
    ],
    junctions: [at.start, at.start + ins.seq.length].map((p, i) => ({
      position: p, type: kind === 'TA' ? 'ta' : 'blunt', overhang: kind === 'TA' ? 'A/T' : '', enzymes: [], leftInput: i === 0 ? 0 : 1, rightInput: i === 0 ? 1 : 0,
      sequence: { top: doc.seq.slice(Math.max(0, p - 6), p + 6), bottom: '', start: Math.max(0, p - 6) }, regeneratedSites: [],
    })),
  };
  addHistory(doc, 'TOPO cloning', `${insert.name} (${orientation}) into ${vector.name} at ${at.start + 1}–${at.end}`);
  const warnings = [];
  if (kind !== 'directional') warnings.push('This kind of cloning is not directional: pick colonies and check the orientation.');
  if (kind === 'TA') warnings.push('The insert needs a 3′ A overhang (Taq-amplified); proofreading polymerases leave blunt ends.');
  return { doc, orientation, warnings };
}
