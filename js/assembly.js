// Assembly beyond two restriction fragments.
//
// * homologyAssemble()  — N fragments joined by shared terminal overlaps
//   (scarless; circular or linear).
// * goldenGateAssemble() — Type IIS fragments with explicit 4-nt overhangs, or
//   fragments whose overhangs are derived from the outermost cut sites.
// * Plans are self-contained JSON (`bublagene/assembly/1`) and replay exactly.
// * batchAssembly() builds many constructs from one named fragment library.
//
// Overhang convention (matches cloning.js): an overhang is written as the
// top-strand bases between the two cuts. Two ends ligate when those strings are
// equal. A fragment's body is the double-stranded part between the inner cuts;
// each overhang appears exactly once in the product.

import { revcomp, subseq, mod } from './seq.js';
import { newDocument, newId, addHistory } from './model.js';
import { gibsonAssemble, CLONING_SCHEMA } from './cloning.js';
import { ENZYME_BY_NAME, findSites, cutEnd } from './enzymes.js';
import { designPrimerPairs } from './primers.js';

export const ASSEMBLY_SCHEMA = 'bublagene/assembly/1';

// Deterministic 64-bit FNV-1a (two 32-bit lanes) for plan replay. Not a
// cryptographic digest; it exists so a plan can prove it rebuilt the same bases.
export function contentHash(seq) {
  let h1 = 0x811c9dc5;
  let h2 = 0x01000193;
  for (let i = 0; i < seq.length; i++) {
    const c = seq.charCodeAt(i) & 0xff;
    h1 = Math.imul(h1 ^ c, 0x01000193) >>> 0;
    h2 = Math.imul(h2 ^ (c + i), 0x85ebca6b) >>> 0;
  }
  return `${h1.toString(16).padStart(8, '0')}${h2.toString(16).padStart(8, '0')}`;
}

function baseName(name) {
  return String(name ?? 'fragment').replace(/\.(dna|gb|gbk|genbank|fasta|fa|fna)$/i, '');
}

function inputDoc(input) {
  if (typeof input === 'string') return { name: 'fragment', seq: input.toUpperCase(), circular: false, features: [], primers: [] };
  if (input?.doc) return input.doc;
  if (input?.seq) return input;
  throw new Error('Each fragment needs a sequence (a string, a document or { seq })');
}

// ---- homology assembly ------------------------------------------------------

export function homologyAssemble(inputs, {
  circular = true, minOverlap = 15, name = null, description = '', maxFragments = 12,
} = {}) {
  const result = gibsonAssemble(inputs, { circular, minOverlap, name, description, maxFragments });
  const plan = planFromGibson(result.doc, inputs, { circular, minOverlap });
  result.doc.assembly = plan;
  result.plan = plan;
  return result;
}

function planFromGibson(doc, inputs, { circular, minOverlap }) {
  const c = doc.cloning;
  const fragments = c.inputs.map((i, index) => {
    const src = inputDoc(inputs[i.index]);
    return {
      index,
      name: baseName(i.name),
      role: i.role,
      orientation: i.orientation,
      seqLen: src.seq.length,
      digest: contentHash(src.seq),
      seq: src.seq,
      left: null,
      right: null,
      overlapIn: i.overlapIn ?? 0,
      overlapOut: i.overlapOut ?? 0,
      retained: i.fragment,
    };
  });
  return {
    schema: ASSEMBLY_SCHEMA,
    method: 'homology',
    created: new Date().toISOString(),
    circular,
    minOverlap,
    enzyme: null,
    product: { name: doc.name, seqLen: doc.seq.length, digest: contentHash(doc.seq) },
    fragments,
    junctions: c.junctions.map((j, n) => ({
      position: j.position,
      type: 'overlap',
      overhang: j.overhang,
      scarLength: 0,
      leftFragment: j.leftInput,
      rightFragment: j.rightInput,
      evidence: { sequence: j.sequence, regeneratedSites: j.regeneratedSites, index: n },
    })),
    history: doc.history?.map((h) => ({ when: h.when, action: h.action, detail: h.detail })) ?? [],
  };
}

// ---- golden gate ------------------------------------------------------------

// Derive {body, left, right} from the outermost Type IIS cuts of a document.
// The left site must point right (strand +1) and the right site must point left
// (strand -1). With `circular`, the sequence is treated as one period.
export function goldenGatePieces(input, { enzyme = 'BsaI', leftSite = null, rightSite = null } = {}) {
  const doc = inputDoc(input);
  const e = typeof enzyme === 'string' ? ENZYME_BY_NAME[enzyme] : enzyme;
  if (!e) throw new Error(`Unknown enzyme ${enzyme}`);
  const L = doc.seq.length;
  const sites = findSites(doc.seq, !!doc.circular, e);
  const forward = sites.filter((s) => s.strand === 1);
  const reverse = sites.filter((s) => s.strand === -1);
  const left = leftSite != null ? (typeof leftSite === 'number' ? sites[leftSite] : leftSite) : forward.sort((a, b) => a.top - b.top)[0];
  const right = rightSite != null ? (typeof rightSite === 'number' ? sites[rightSite] : rightSite) : reverse.sort((a, b) => b.top - a.top)[0];
  if (!left || !right) throw new Error(`${doc.name}: needs one forward and one reverse ${e.name} site`);
  const leftEnd = cutEnd(doc.seq, !!doc.circular, left);
  const rightEnd = cutEnd(doc.seq, !!doc.circular, right);
  const leftStart = Math.min(left.top, left.bottom);
  const leftBody = Math.max(left.top, left.bottom);
  const rightBody = Math.min(right.top, right.bottom);
  const rightStop = Math.max(right.top, right.bottom);
  const body = doc.circular
    ? subseq(doc.seq, leftBody, rightBody, true)
    : doc.seq.slice(leftBody, rightBody);
  if (!body.length) throw new Error(`${doc.name}: the two cut sites leave no body`);
  return {
    name: baseName(doc.name),
    seq: body.toUpperCase(),
    left: leftEnd.oh,
    right: rightEnd.oh,
    doc,
    sites: { left, right },
    evidence: {
      left: { pos: left.top, overhang: leftEnd.oh, type: leftEnd.type },
      right: { pos: right.top, overhang: rightEnd.oh, type: rightEnd.type },
      bodyRange: { start: leftBody, end: rightBody },
      outerRange: { start: leftStart, end: rightStop },
    },
  };
}

function explicitFragment(input, index) {
  if (input?.left != null && input?.right != null && input?.seq) {
    return {
      name: baseName(input.name ?? `fragment ${index + 1}`),
      seq: String(input.seq).toUpperCase(),
      left: String(input.left).toUpperCase(),
      right: String(input.right).toUpperCase(),
      doc: input.doc ?? null,
      orientation: input.orientation ?? 'forward',
      evidence: input.evidence ?? null,
    };
  }
  if (input?.doc || typeof input === 'string' || input?.seq) {
    return goldenGatePieces(input, { enzyme: input.enzyme });
  }
  throw new Error('Golden Gate fragments need { seq, left, right } or a document with flanking Type IIS sites');
}

function fragmentFromInput(input, index, opts) {
  if (input?.left != null && input?.right != null) return explicitFragment(input, index);
  if (input?.doc || typeof input === 'string' || input?.seq) {
    const pieces = goldenGatePieces(input, { enzyme: opts.enzyme });
    return pieces;
  }
  return explicitFragment(input, index);
}

// Join Type IIS fragments. Each fragment contributes left overhang + body; the
// next fragment's left overhang is the shared junction. For a circular product
// the last fragment's right overhang must equal the first fragment's left.
export function goldenGateAssemble(fragments, {
  circular = true, enzyme = 'BsaI', name = null, description = '',
  scar = 'minimal', forbiddenSites = true, allowDuplicateOverhangs = false, maxFragments = 16,
} = {}) {
  if (!Array.isArray(fragments) || fragments.length < 2) throw new Error('Golden Gate needs at least two fragments');
  if (fragments.length > maxFragments) throw new Error(`Golden Gate supports at most ${maxFragments} fragments here`);
  const pieces = fragments.map((f, i) => fragmentFromInput(f, i, { enzyme }));
  for (const [i, p] of pieces.entries()) {
    if (!p.left || !p.right) throw new Error(`Fragment ${p.name} needs both overhangs`);
    if (p.left.length !== p.right.length) throw new Error(`Fragment ${p.name} overhangs differ in length`);
  }
  const overhangLength = pieces[0].left.length;
  if (!overhangLength) throw new Error('Golden Gate overhangs cannot be empty');
  for (const p of pieces) {
    if (p.left.length !== overhangLength || p.right.length !== overhangLength) {
      throw new Error(`Fragment ${p.name} overhangs must all be ${overhangLength} nt`);
    }
  }
  const junctions = [];
  for (let i = 0; i < pieces.length; i++) {
    const a = pieces[i];
    const b = pieces[(i + 1) % pieces.length];
    if (!circular && i === pieces.length - 1) break;
    if (a.right !== b.left) {
      throw new Error(`Overhangs do not match: ${a.name} ends ${a.right}, ${b.name} starts ${b.left}`);
    }
    junctions.push({
      left: a.right, right: b.left, compatible: true,
      leftFragment: i, rightFragment: (i + 1) % pieces.length,
      scar: a.right, scarLength: a.right.length,
    });
  }
  if (circular && pieces[pieces.length - 1].right !== pieces[0].left) {
    throw new Error(`The circle does not close: ${pieces[pieces.length - 1].name} ends ${pieces[pieces.length - 1].right}, ${pieces[0].name} starts ${pieces[0].left}`);
  }

  let seq = '';
  const placements = [];
  for (let i = 0; i < pieces.length; i++) {
    const p = pieces[i];
    const start = seq.length;
    seq += p.left + p.seq;
    placements.push({ index: i, start, end: seq.length, piece: p });
  }
  if (!circular) seq += pieces[pieces.length - 1].right;
  // For a circle the first overhang closes the last right overhang and is
  // already present at the start, so nothing is appended.

  const productName = name ?? `${pieces[0].name} Golden Gate`;
  const doc = newDocument({ name: productName, seq, circular });
  doc.features = [];
  for (const placement of placements) {
    if (!placement.piece.doc) continue;
    const src = placement.piece.doc;
    const offset = placement.start + placement.piece.left.length;
    for (const f of src.features ?? []) {
      const segs = f.segments
        .map((s) => ({ ...structuredClone(s), start: s.start - 0 + offset, end: s.end - 0 + offset }))
        .filter((s) => s.start < placement.end && s.end > placement.start);
      if (!segs.length) continue;
      doc.features.push({ ...structuredClone(f), id: newId(), segments: segs });
    }
  }

  const plan = {
    schema: ASSEMBLY_SCHEMA,
    method: 'golden-gate',
    created: new Date().toISOString(),
    circular,
    enzyme: typeof enzyme === 'string' ? enzyme : /** @type {any} */ (enzyme).name,
    scar,
    product: { name: doc.name, seqLen: doc.seq.length, digest: contentHash(doc.seq) },
    fragments: pieces.map((p, index) => ({
      index, name: p.name, role: index === 0 ? 'backbone' : 'insert',
      orientation: p.orientation ?? 'forward', seqLen: p.seq.length + p.left.length + p.right.length,
      digest: contentHash(p.seq), seq: p.seq, left: p.left, right: p.right,
      evidence: p.evidence ?? null,
    })),
    junctions: junctions.map((j, n) => {
      const at = placements[j.rightFragment]?.start ?? 0;
      const window = circular ? subseq(seq, mod(at - 6, seq.length || 1), at + 6, true) : seq.slice(Math.max(0, at - 6), at + 6);
      return {
        position: at, type: 'overhang', overhang: j.left, scar: j.scar, scarLength: j.scarLength,
        leftFragment: j.leftFragment, rightFragment: j.rightFragment,
        evidence: { index: n, window, regeneratedSites: [] },
      };
    }),
    history: [],
  };

  const checks = validateJunctions({ doc, plan, enzyme, forbiddenSites, allowDuplicateOverhangs, frame: false });
  plan.checks = checks;
  doc.assembly = plan;
  doc.cloning = {
    schema: CLONING_SCHEMA, type: 'golden-gate', operation: 'golden-gate-assembly',
    created: plan.created, circular,
    product: { name: doc.name, seqLen: doc.seq.length },
    inputs: plan.fragments.map((f) => ({ name: f.name, seqLen: f.seqLen, role: f.role, orientation: f.orientation, ends: { left: { type: '5', oh: f.left }, right: { type: '5', oh: f.right } } })),
    junctions: plan.junctions.map((j) => ({ position: j.position, type: 'overhang', overhang: j.overhang, enzymes: [], leftInput: j.leftFragment, rightInput: j.rightFragment, sequence: j.evidence.window, regeneratedSites: j.evidence.regeneratedSites })),
  };
  addHistory(doc, 'Golden Gate assembly', `${pieces.map((p) => p.name).join(' + ')} → ${doc.seq.length} bp ${circular ? 'circular' : 'linear'}; ${junctions.length} junction(s)`);
  doc.notes.description = description || `Golden Gate assembly with ${plan.enzyme}`;
  const errors = checks.filter((c) => c.severity === 'error');
  return { doc, plan, junctions: plan.junctions, checks, warnings: checks.filter((c) => c.severity === 'warning').map((c) => c.message), errors };
}

// ---- junction validation ----------------------------------------------------

/** @param {{doc?:any,plan?:any,enzyme?:any,forbiddenSites?:boolean,allowDuplicateOverhangs?:boolean,frame?:any}} [options] */
export function validateJunctions({ doc, plan, enzyme = null, forbiddenSites = true, allowDuplicateOverhangs = false, frame = false } = {}) {
  const checks = [];
  const seen = new Map();
  for (const [n, j] of plan.junctions.entries()) {
    if (j.overhang) {
      if (seen.has(j.overhang) && !allowDuplicateOverhangs) {
        checks.push({ severity: 'warning', code: 'duplicate-overhang', junction: n, message: `Overhang ${j.overhang} is used at junctions ${seen.get(j.overhang) + 1} and ${n + 1}; misassembly is possible.` });
      }
      if (!seen.has(j.overhang)) seen.set(j.overhang, n);
    }
  }
  const names = plan.fragments.map((f) => f.name);
  const dupes = names.filter((n, i) => names.indexOf(n) !== i);
  if (dupes.length) checks.push({ severity: 'warning', code: 'duplicate-fragment', message: `Duplicate fragment name(s): ${[...new Set(dupes)].join(', ')}` });
  if (forbiddenSites && enzyme && ENZYME_BY_NAME[enzyme]) {
    const sites = findSites(doc.seq, doc.circular, ENZYME_BY_NAME[enzyme]);
    if (sites.length) {
      checks.push({
        severity: 'error', code: 'forbidden-site',
        message: `${sites.length} internal ${enzyme} site(s) remain; the assembly would be re-cut.`,
        positions: sites.map((s) => s.top),
      });
    }
  }
  if (frame) {
    const frameCheck = checkReadingFrame(doc, frame);
    if (frameCheck) checks.push(frameCheck);
  }
  return checks;
}

// Check that every CDS crossing a junction keeps a continuous reading frame.
export function checkReadingFrame(doc, junctions) {
  const positions = (junctions ?? doc.assembly?.junctions ?? []).map((j) => j.position).sort((a, b) => a - b);
  for (const f of doc.features ?? []) {
    if (f.type !== 'CDS') continue;
    const span = f.segments[0];
    const crosses = positions.some((p) => p > span.start && p < span.end);
    if (!crosses) continue;
    const len = f.segments.reduce((n, s) => n + Math.max(0, s.end - s.start), 0);
    if (len % 3 !== 0) {
      return { severity: 'error', code: 'frame-shift', feature: f.name, message: `${f.name} crosses a junction and is not a whole number of codons (${len} bp).` };
    }
  }
  return null;
}

// ---- plans ------------------------------------------------------------------

export function assemblyPlan(doc) {
  const p = doc?.assembly;
  if (!p) return null;
  return structuredClone(p);
}

export function assemblyFromJSON(value) {
  const plan = typeof value === 'string' ? JSON.parse(value) : structuredClone(value);
  if (!plan || plan.schema !== ASSEMBLY_SCHEMA) throw new Error(`Unsupported assembly plan (expected schema ${ASSEMBLY_SCHEMA})`);
  if (!plan.product || !Number.isFinite(plan.product.seqLen)) throw new Error('Assembly plan is missing the product size');
  if (!Array.isArray(plan.fragments) || !Array.isArray(plan.junctions)) throw new Error('Assembly plan is missing fragments or junctions');
  return plan;
}

// Rebuild the product from the plan alone and confirm the digest.
export function replayAssembly(plan) {
  const p = assemblyFromJSON(plan);
  let seq;
  if (p.method === 'golden-gate') {
    seq = p.fragments.map((f) => (f.left ?? '') + f.seq).join('');
    if (!p.circular) seq += p.fragments[p.fragments.length - 1].right ?? '';
  } else {
    // homology plans store retained regions; replay via overlap walk
    seq = p.fragments[0].seq ?? '';
    for (let i = 1; i < p.fragments.length; i++) {
      const ov = p.fragments[i].overlapIn ?? p.junctions.find((j) => j.rightFragment === i)?.overhang?.length ?? p.minOverlap;
      seq += (p.fragments[i].seq ?? '').slice(ov);
    }
    if (p.circular) seq = seq.slice(0, seq.length - (p.junctions.find((j) => j.position === p.product.seqLen)?.overhang?.length ?? 0));
  }
  const digest = contentHash(seq);
  return { ok: digest === p.product.digest && seq.length === p.product.seqLen, seq, digest, expected: p.product };
}

// ---- assembly primers -------------------------------------------------------

// Design PCR primers that add the overlaps each fragment needs. For linear
// fragments the forward primer carries the previous fragment's 3' end and the
// reverse primer carries the next fragment's 5' end.
export function designAssemblyPrimers(fragments, {
  overlap = 25, targetTm = 60, minLen = 18, maxLen = 32, circular = true, maxTmDiff = 10,
} = {}) {
  const docs = fragments.map(inputDoc);
  return docs.map((doc, i) => {
    const prev = docs[(i - 1 + docs.length) % docs.length];
    const next = docs[(i + 1) % docs.length];
    const fwdTail = circular || i > 0 ? prev.seq.slice(-overlap) : '';
    const revTail = circular || i < docs.length - 1 ? next.seq.slice(0, overlap) : '';
    const design = designPrimerPairs(doc.seq, false, 0, doc.seq.length, {
      targetTm, minLen, maxLen, fwdTail, revTail, maxTmDiff,
    });
    return {
      name: doc.name,
      forward: { sequence: design.forward.seq, tail: fwdTail, tm: design.forward.tm },
      reverse: { sequence: design.reverse.seq, tail: revTail, tm: design.reverse.tm },
      productLength: doc.seq.length + fwdTail.length + revTail.length,
    };
  });
}

// ---- batch assembly ---------------------------------------------------------

// library: { name: doc | seq | { seq, left, right } }
// designs: [{ name, fragments: [name | { name, reverse }], method, ... }]
export function batchAssembly(library, designs, opts = {}) {
  const results = [];
  for (const design of designs) {
    try {
      const inputs = design.fragments.map((f) => {
        const key = typeof f === 'string' ? f : f.name;
        const src = library[key];
        if (!src) throw new Error(`Unknown fragment "${key}"`);
        const item = typeof src === 'string' ? { name: key, seq: src } : { ...src, name: src.name ?? key };
        return typeof f === 'object' && f.reverse && item.doc ? { ...item, doc: { ...item.doc, seq: revcomp(item.doc.seq) } } : item;
      });
      const method = design.method ?? 'homology';
      const result = method === 'golden-gate'
        ? goldenGateAssemble(inputs, { ...opts, ...design, name: design.name })
        : homologyAssemble(inputs, { ...opts, ...design, name: design.name });
      results.push({ name: design.name, ok: true, doc: result.doc, plan: result.plan, warnings: result.warnings ?? [] });
    } catch (e) {
      results.push({ name: design.name, ok: false, error: e.message });
    }
  }
  return { results, ok: results.filter((r) => r.ok).length, failed: results.filter((r) => !r.ok).length };
}
