// Restriction cloning and Gibson assembly.
//
// A digested fragment (see enzymes.js digest()) becomes a *piece*: the double-
// stranded molecule that comes out of the digest. The piece's top strand is what
// a product inherits when the piece goes in forward; when the piece is flipped,
// the product inherits its bottom strand (reverse complemented). That distinction
// matters: an EcoRI piece with a 5' overhang has a top strand four bases longer
// than its bottom strand, so a flipped piece contributes the shorter strand.
//
// Ends are { type, oh, enzyme }. `oh` is the overhang written as top-strand
// bases, so EcoRI (G^AATT_C) is { type: '5', oh: 'AATT' } and BamHI
// (G^GATC_C) is { type: '5', oh: 'GATC' }. Two ends ligate when their types
// match and their overhangs are equal (EcoRI with EcoRI, BamHI with BglII,
// blunt with blunt). Natural (uncut) ends do not ligate.
//
// Coordinates are 0-based end-exclusive everywhere, as in the rest of BublaGene.
//
// A construct keeps `doc.cloning`, a structured record of what it was built
// from (operation, parents, fragments and junctions). `historyTreeXML()` renders
// it in the binary format's history-tree schema.

import { revcomp, mod, subseq, rangeLength, inRange } from './seq.js';
import { newDocument, newId, addHistory } from './model.js';
import { findSites, ENZYME_BY_NAME } from './enzymes.js';
import { el, toXML } from './xml.js';

const FLANK = 12;

// ---- ends ------------------------------------------------------------------

export function describeEnd(end) {
  if (!end) return 'unknown end';
  if (end.type === 'blunt') return 'blunt end';
  if (end.type === 'natural') return 'uncut end';
  if (end.type === 'overlap') return `${end.oh.length} bp overlap`;
  return `${end.type}′ overhang ${end.oh}`;
}

// Can these two ends be ligated? Equal top-strand overhang strings, not reverse
// complements: one side's overhang is on the bottom strand.
export function endsCompatible(a, b) {
  if (!a || !b) return false;
  if (a.type === 'blunt' && b.type === 'blunt') return true;
  if (a.type !== b.type) return false;
  if (a.type !== '5' && a.type !== '3') return false;
  if (!a.oh || !b.oh || a.oh.length !== b.oh.length) return false;
  return a.oh.toUpperCase() === b.oh.toUpperCase();
}

const endOf = ({ cut, end }) => ({ type: end.type, oh: end.oh ?? '', enzyme: cut?.enzyme?.name ?? null });

const flipEnd = (end, enzyme) =>
  end.type === 'natural'
    ? { type: 'natural', oh: '', enzyme: null }
    : { type: end.type, oh: revcomp(end.oh), enzyme: enzyme ?? null };

function span(a, b, L, circular) {
  return circular ? mod(b - a, L) || L : b - a;
}

// ---- pieces ----------------------------------------------------------------

// Which parts of a source segment [start, start+sl) fall inside the piece
// region [from, from+len)? The source may be circular, so a segment can cross
// the region boundary. Returns intervals relative to the piece.
function intersectRegion(start, sl, from, len, L, circular) {
  if (!circular) {
    const a = Math.max(0, start - from);
    const b = Math.min(len, start + sl - from);
    return b > a ? [[a, b]] : [];
  }
  const r = mod(start - from, L);
  const pieces = [];
  if (r < len) pieces.push([r, Math.min(r + sl, len)]);
  if (r + sl > L) {
    const tail = r + sl - L;
    if (tail > 0) {
      const b = Math.min(tail, len);
      if (b > 0) pieces.push([0, b]);
    }
  }
  return pieces;
}

function mapFeatures(doc, from, len, circular, reverse) {
  const L = doc.seq.length;
  const out = [];
  for (const f of doc.features ?? []) {
    const pieces = [];
    (f.segments ?? []).forEach((s, srcSeg) => {
      const sl = rangeLength(s, L);
      for (const [a, b] of intersectRegion(s.start, sl, from, len, L, circular)) {
        pieces.push({
          seg: s,
          srcSeg,
          start: reverse ? len - b : a,
          end: reverse ? len - a : b,
        });
      }
    });
    if (!pieces.length) continue;
    pieces.sort((x, y) => x.start - y.start);
    // Segments of one feature stay together across designed gaps (introns), but
    // a single segment cut by the enzyme becomes separate annotations.
    const runs = [];
    for (const p of pieces) {
      const run = runs[runs.length - 1];
      const tail = run?.[run.length - 1];
      if (run && (tail.end === p.start || !run.some((q) => q.srcSeg === p.srcSeg))) run.push(p);
      else runs.push([p]);
    }
    for (const run of runs) {
      out.push({
        ...structuredClone(f),
        id: newId(),
        strand: reverse ? -f.strand : f.strand,
        segments: run.map(({ seg, start, end }) => ({ ...structuredClone(seg), start, end })),
      });
    }
  }
  return out;
}

// Build the piece a digested fragment represents. `reverse` flips it, taking the
// bottom strand as the piece's top strand and reversing its features.
export function pieceFromFragment(doc, fragment, { reverse = false, name = null } = {}) {
  if (!doc?.seq || !fragment?.left?.cut || !fragment?.right?.cut) {
    throw new Error('pieceFromFragment needs a document and a fragment from digest()');
  }
  const L = doc.seq.length;
  const circular = !!doc.circular;
  const from = reverse ? fragment.left.cut.bottom : fragment.left.cut.top;
  const to = reverse ? fragment.right.cut.bottom : fragment.right.cut.top;
  const len = span(from, to, L, circular);
  const raw = subseq(doc.seq, from, to, circular);
  const seq = (reverse ? revcomp(raw) : raw).toUpperCase();
  const features = mapFeatures(doc, from, len, circular, reverse);
  const left = reverse ? flipEnd(endOf(fragment.right), fragment.right.cut.enzyme?.name) : endOf(fragment.left);
  const right = reverse ? flipEnd(endOf(fragment.left), fragment.left.cut.enzyme?.name) : endOf(fragment.right);
  return {
    doc: { name: name ?? doc.name, seq, circular: false, features, primers: [] },
    source: doc,
    fragment,
    reverse,
    left,
    right,
  };
}

function asPiece(x, reverse = false) {
  if (!x) throw new Error('Missing fragment');
  if (x.source && x.left && x.right && x.doc) {
    if (!reverse || x.reverse === reverse) return x;
    if (!x.fragment) throw new Error('Cannot flip this piece: it has no source fragment');
    return pieceFromFragment(x.source, x.fragment, { reverse, name: x.doc.name });
  }
  if (x.doc && x.fragment) return pieceFromFragment(x.doc, x.fragment, { reverse });
  throw new Error('Expected { doc, fragment } or a piece from pieceFromFragment()');
}

function junctionBetween(leftPiece, rightPiece, position) {
  const a = leftPiece.right;
  const b = rightPiece.left;
  const ok = endsCompatible(a, b);
  const names = [a.enzyme, b.enzyme].filter(Boolean);
  const problem = ok
    ? null
    : `${describeEnd(a)} cannot join ${describeEnd(b)}${names.length ? ` (${names.join('/')})` : ''}`;
  return { position, left: a, right: b, compatible: ok, type: a.type, oh: a.oh, enzymes: names, problem };
}

// ---- checking --------------------------------------------------------------

// Test a vector backbone and an insert in both orientations, without building
// anything. Answers "will it ligate?", "can it go in backwards?", and sizes.
export function checkLigation({ backbone, insert, circular = true } = {}) {
  const b = asPiece(backbone, false);
  const variants = {};
  for (const orientation of ['forward', 'reverse']) {
    const ins = asPiece(insert, orientation === 'reverse');
    const junctions = circular
      ? [junctionBetween(b, ins, b.doc.seq.length), junctionBetween(ins, b, b.doc.seq.length + ins.doc.seq.length)]
      : [junctionBetween(b, ins, b.doc.seq.length)];
    const problems = junctions.filter((j) => !j.compatible).map((j) => j.problem);
    variants[orientation] = {
      orientation,
      compatible: problems.length === 0,
      junctions,
      problems,
      productSize: b.doc.seq.length + ins.doc.seq.length,
    };
  }
  const orientations = ['forward', 'reverse'].filter((o) => variants[o].compatible);
  const warnings = [];
  if (orientations.length === 2) {
    warnings.push('The insert ligates in either orientation; pick one if direction matters.');
  }
  return {
    circular,
    forward: variants.forward,
    reverse: variants.reverse,
    orientations,
    preferred: orientations.includes('forward') ? 'forward' : orientations[0] ?? null,
    warnings,
  };
}

// ---- product ---------------------------------------------------------------

function junctionContext(seq, circular, pos, flank) {
  const L = seq.length;
  const want = Math.max(1, Math.min(2 * flank, L));
  if (circular && want === L) {
    const start = mod(pos, L);
    return { start, top: subseq(seq, start, start, true) };
  }
  if (!circular) {
    const start = Math.max(0, Math.min(L - want, pos - flank));
    return { start, top: seq.slice(start, start + want) };
  }
  const start = mod(pos - flank, L);
  return { start, top: subseq(seq, start, start + want, true) };
}

function regeneratedSites(seq, circular, pos, names) {
  const L = seq.length;
  const out = [];
  const seen = new Set();
  for (const name of [...new Set((names ?? []).filter(Boolean))]) {
    const enzyme = ENZYME_BY_NAME[name];
    if (!enzyme) continue;
    for (const site of findSites(seq, circular, enzyme)) {
      if (!inRange(mod(pos, L), site.site, L)) continue;
      const key = `${name}:${site.top}`;
      if (seen.has(key)) continue;
      seen.add(key);
      out.push({ name, pos: site.top, cut: site.top, start: site.site.start, end: site.site.end, site: enzyme.site, siteCount: 1 });
    }
  }
  return out.sort((a, b) => a.pos - b.pos);
}

function baseName(name) {
  return String(name ?? 'fragment').replace(/\.(dna|gb|gbk|genbank|fasta|fa|fna)$/i, '');
}

function defaultName(pieces, type) {
  const first = baseName(pieces[0]?.source?.name ?? pieces[0]?.doc?.name);
  if (pieces.length === 1) return first;
  if (type === 'gibson') return `${first} Gibson`;
  const inserts = pieces.slice(1).map((p) => baseName(p.source?.name ?? p.doc.name));
  return `${first} + ${inserts.join(' + ')}`;
}

// A short provenance line for Notes.description, which (unlike the HistoryTree)
// is plain XML and survives a .dna round trip.
function provenanceText(inputs, type, circular) {
  const label = type === 'gibson' ? 'Gibson assembly' : 'Restriction cloning';
  const parts = inputs.map((i) => {
    const enzymes = [...new Set([i.ends?.left?.enzyme, i.ends?.right?.enzyme].filter(Boolean))];
    const site = enzymes.length ? ` (${enzymes.join('/')})` : '';
    return `${i.name ?? 'fragment'}${site}${i.orientation === 'reverse' ? ' reversed' : ''}`;
  });
  return `${label}: ${parts.join(' + ')}${circular ? '' : ' (linear)'}`;
}

// Join pieces in the given order. Every junction must be compatible. Used by
// restriction cloning; Gibson assembly builds its product in gibsonAssemble().
export function assemblePieces(pieces, {
  circular = false, name = null, type = 'restriction', operation = null, description = '',
} = {}) {
  if (!Array.isArray(pieces) || !pieces.length) throw new Error('Nothing to assemble');
  const n = pieces.length;
  const junctions = [];
  const total = pieces.reduce((sum, p) => sum + p.doc.seq.length, 0);
  let position = 0;
  for (let i = 0; i < n; i++) {
    const j = i < n - 1 || circular ? i : -1;
    if (j < 0) break;
    const next = pieces[(i + 1) % n];
    const at = position + pieces[i].doc.seq.length;
    const junction = junctionBetween(pieces[i], next, at);
    junctions.push({ ...junction, leftInput: i, rightInput: (i + 1) % n });
    position += pieces[i].doc.seq.length;
  }
  const bad = junctions.filter((j) => !j.compatible);
  if (bad.length) throw new Error(`Cannot join the fragments: ${bad.map((j) => j.problem).join('; ')}`);

  const seq = pieces.map((p) => p.doc.seq).join('');
  const doc = newDocument({ name: name ?? defaultName(pieces, type), seq, circular });
  doc.features = [];

  let offset = 0;
  const inputs = [];
  for (let i = 0; i < n; i++) {
    const p = pieces[i];
    const len = p.doc.seq.length;
    for (const f of p.doc.features ?? []) {
      doc.features.push({
        ...structuredClone(f),
        id: newId(),
        segments: f.segments.map((s) => ({ ...structuredClone(s), start: s.start + offset, end: s.end + offset })),
      });
    }
    inputs.push({
      name: p.source?.name ?? p.doc.name,
      seqLen: p.source?.seq.length ?? len,
      circular: !!p.source?.circular,
      role: p.role ?? (i === 0 ? 'backbone' : 'insert'),
      orientation: p.reverse ? 'reverse' : 'forward',
      fragment: p.fragment ? { start: p.fragment.start, end: p.fragment.end } : null,
      piece: { start: offset, end: offset + len, len },
      ends: { left: p.left, right: p.right },
      features: structuredClone(p.source?.features ?? p.doc.features ?? []),
    });
    offset += len;
  }

  for (const j of junctions) {
    const ctx = junctionContext(seq, circular, j.position, FLANK);
    j.sequence = { top: ctx.top, bottom: revcomp(ctx.top), start: ctx.start };
    j.regeneratedSites = regeneratedSites(seq, circular, j.position, j.enzymes);
  }

  const op = operation ?? (type === 'gibson' ? 'gibsonAssembly' : n === 1 && circular ? 'circularize' : 'insertFragment');
  doc.cloning = {
    type,
    operation: op,
    created: new Date().toISOString(),
    circular,
    product: { name: doc.name, seqLen: total },
    inputs,
    junctions: junctions.map((j) => ({
      position: j.position,
      type: j.type,
      overhang: j.oh,
      enzymes: j.enzymes,
      leftInput: j.leftInput,
      rightInput: j.rightInput,
      sequence: j.sequence,
      regeneratedSites: j.regeneratedSites,
    })),
  };

  const what = type === 'gibson' ? 'Gibson assembly' : 'Restriction cloning';
  const labels = inputs.map((i) => `${i.name}${i.orientation === 'reverse' ? ' (reverse)' : ''}`);
  addHistory(
    doc,
    what,
    `${labels.join(' + ')} → ${seq.length} bp ${circular ? 'circular' : 'linear'}; ` +
      `${junctions.length} junction${junctions.length === 1 ? '' : 's'}`,
  );
  doc.notes.description = description || provenanceText(inputs, type, circular);
  return { doc, junctions: doc.cloning.junctions, inputs, warnings: [] };
}

// Convenience wrapper: pick the insert orientation, then ligate. Fragment ends
// must be cut; use pieceFromFragment() or { doc, fragment } from digest().
export function restrictionClone({
  backbone, insert, circular = true, orientation = 'auto', name = null, description = '',
} = {}) {
  const b = { ...asPiece(backbone, false), role: 'backbone' };
  const check = checkLigation({ backbone: b, insert, circular });
  let chosen;
  if (orientation === 'forward' || orientation === 'reverse') {
    if (!check[orientation].compatible) {
      const other = orientation === 'forward' ? 'reverse' : 'forward';
      const extra = check[other].compatible ? ` The insert does ligate in the ${other} orientation.` : '';
      throw new Error(`The insert cannot go in ${orientation}: ${check[orientation].problems.join('; ')}.${extra}`);
    }
    chosen = orientation;
  } else {
    if (!check.orientations.length) {
      throw new Error(`The ends do not ligate: ${check.forward.problems.join('; ')}`);
    }
    chosen = check.preferred;
  }
  const ins = { ...asPiece(insert, chosen === 'reverse'), role: 'insert' };
  const result = assemblePieces([b, ins], {
    circular,
    name: name ?? defaultName([b, ins], 'restriction'),
    type: 'restriction',
    operation: 'insertFragment',
    description,
  });
  result.orientation = chosen;
  result.check = check;
  result.warnings.push(...check.warnings);
  if (check.orientations.length === 1 && chosen === 'reverse') {
    result.warnings.push('The insert only ligates in the reverse orientation.');
  }
  return result;
}

// ---- Gibson assembly -------------------------------------------------------

function gibsonPiece(input, reverse) {
  let src;
  if (typeof input === 'string') {
    src = { name: 'fragment', seq: input, circular: false, features: [], primers: [] };
  } else if (input?.doc && Number.isFinite(input.start) && Number.isFinite(input.end)) {
    const d = input.doc;
    const region = subseq(d.seq, input.start, input.end, d.circular);
    const features = mapFeatures(d, input.start, region.length, !!d.circular, false);
    src = { name: input.name ?? d.name, seq: region, circular: false, features, primers: [] };
  } else if (input?.doc) {
    src = input.doc;
  } else {
    src = input;
  }
  if (!src?.seq || typeof src.seq !== 'string') {
    throw new Error('Each Gibson fragment needs a sequence (a document or { seq })');
  }
  if (src.circular) {
    throw new Error(`"${src.name}" is circular; Gibson needs linear fragments (cut or PCR-amplified first)`);
  }
  const seq = src.seq.toUpperCase();
  const L = seq.length;
  const features = [];
  for (const f of src.features ?? []) {
    const segs = reverse
      ? [...f.segments].reverse().map((s) => ({ ...structuredClone(s), start: L - s.end, end: L - s.start }))
      : f.segments.map((s) => structuredClone(s));
    features.push({ ...structuredClone(f), id: newId(), strand: reverse ? -f.strand : f.strand, segments: segs });
  }
  return {
    doc: { name: src.name ?? 'fragment', seq: reverse ? revcomp(seq) : seq, circular: false, features, primers: [] },
    source: src,
    fragment: null,
    reverse,
    left: { type: 'overlap', oh: '', enzyme: null },
    right: { type: 'overlap', oh: '', enzyme: null },
  };
}

function maxOverlap(a, b, min, cap) {
  const most = Math.min(a.length, b.length, cap ?? Infinity);
  for (let k = most; k >= min; k--) if (a.endsWith(b.slice(0, k))) return k;
  return 0;
}

function clipFeatures(features, start, end) {
  const out = [];
  for (const f of features ?? []) {
    const segs = [];
    for (const s of f.segments ?? []) {
      const a = Math.max(s.start, start);
      const b = Math.min(s.end, end);
      if (b > a) segs.push({ ...structuredClone(s), start: a, end: b });
    }
    if (segs.length) out.push({ ...structuredClone(f), id: newId(), segments: segs });
  }
  return out;
}

function orientationKey(placements, count) {
  const key = new Array(count).fill(0);
  for (const p of placements) key[p.index] = p.piece.reverse ? 1 : 0;
  return key.join('');
}

// Rotations of a circular assembly describe the same molecule, and reversing a
// whole cycle is its reverse-complement strand. Canonicalise both before
// comparing solutions.
function cycleKey(solution, circular) {
  const parts = solution.placements.map((p) => ({ index: p.index, reverse: p.piece.reverse ? 1 : 0 }));
  const linearKey = (arr) => arr.map((p) => `${p.index}:${p.reverse}`).join(',');
  if (!circular || parts.length < 2) return linearKey(parts);
  const canonical = (arr) => {
    let best = 0;
    for (let i = 1; i < arr.length; i++) if (arr[i].index < arr[best].index) best = i;
    return linearKey([...arr.slice(best), ...arr.slice(0, best)]);
  };
  const forward = canonical(parts);
  const reverse = canonical([...parts].reverse().map((p) => ({ index: p.index, reverse: p.reverse ? 0 : 1 })));
  return forward < reverse ? forward : reverse;
}

const totalOverlap = (solution) =>
  solution.placements.reduce((sum, p) => sum + p.overlapIn, 0) + solution.closing;

const reversals = (solution) => solution.placements.filter((p) => p.piece.reverse).length;

// Assemble linear fragments that share terminal overlaps. The standard is 15–40
// bp of identical sequence: each fragment's 3' end overlaps the next one's 5'
// end. Order and orientation are inferred unless the inputs are already in order.
export function gibsonAssemble(inputs, {
  circular = true, minOverlap = 15, name = null, description = '', maxFragments = 8,
} = {}) {
  if (!Array.isArray(inputs) || !inputs.length) throw new Error('Nothing to assemble');
  if (inputs.length > maxFragments) {
    throw new Error(`Gibson assembly supports at most ${maxFragments} fragments here (${inputs.length} given)`);
  }
  if (!(minOverlap >= 1)) throw new Error('minOverlap must be at least 1');
  if (!circular && inputs.length < 2) throw new Error('A linear Gibson assembly needs at least two fragments');

  const candidates = inputs.map((input, index) => ({
    index,
    forward: gibsonPiece(input, false),
    reverse: gibsonPiece(input, true),
  }));

  const finished = [];
  let visited = 0;
  const LIMIT = 50000;

  function dfs(acc, placements, used) {
    if (++visited > LIMIT) return;
    if (used.size === candidates.length) {
      let closing = 0;
      if (circular) {
        closing = maxOverlap(acc, placements[0].piece.doc.seq, minOverlap, acc.length - 1);
        if (!closing) return;
      }
      finished.push({ acc, placements: placements.map((p) => ({ ...p })), closing });
      return;
    }
    for (let j = 0; j < candidates.length; j++) {
      if (used.has(j)) continue;
      for (const orientation of ['forward', 'reverse']) {
        const piece = candidates[j][orientation];
        const ov = maxOverlap(acc, piece.doc.seq, minOverlap, piece.doc.seq.length - 1);
        if (!ov) continue;
        const offset = acc.length - ov;
        used.add(j);
        placements.push({ piece, index: j, overlapIn: ov, offset });
        dfs(acc + piece.doc.seq.slice(ov), placements, used);
        placements.pop();
        used.delete(j);
      }
    }
  }

  for (let start = 0; start < candidates.length; start++) {
    if (inputs.length === 1) {
      // A single fragment can only circularize on itself.
      const piece = candidates[start].forward;
      const closing = maxOverlap(piece.doc.seq, piece.doc.seq, minOverlap, piece.doc.seq.length - 1);
      if (circular && closing) {
        finished.push({ acc: piece.doc.seq, placements: [{ piece, index: start, overlapIn: 0, offset: 0 }], closing });
      }
      continue;
    }
    for (const orientation of ['forward', 'reverse']) {
      const piece = candidates[start][orientation];
      const used = new Set([start]);
      dfs(piece.doc.seq, [{ piece, index: start, overlapIn: 0, offset: 0 }], used);
    }
  }

  if (!finished.length) {
    throw new Error(
      `No assembly found: no ends share at least ${minOverlap} matching bases in the same order. ` +
        `Check the fragment orientation or lower minOverlap.`,
    );
  }
  finished.sort(
    (a, b) =>
      totalOverlap(b) - totalOverlap(a) ||
      reversals(a) - reversals(b) ||
      orientationKey(a.placements, candidates.length).localeCompare(orientationKey(b.placements, candidates.length)) ||
      a.placements[0].index - b.placements[0].index,
  );
  const best = finished[0];
  const placements = best.placements;
  const productSeq = circular ? best.acc.slice(0, best.acc.length - best.closing) : best.acc;
  const productLen = productSeq.length;
  const n = placements.length;

  const doc = newDocument({ name: name ?? defaultName(placements.map((p) => p.piece), 'gibson'), seq: productSeq, circular });
  doc.features = [];
  const inputsMeta = [];
  const junctions = [];

  for (let i = 0; i < n; i++) {
    const { piece, overlapIn, offset, index } = placements[i];
    const N = piece.doc.seq.length;
    const nextOverlap = i < n - 1 ? placements[i + 1].overlapIn : circular ? best.closing : 0;
    const retainStart = i === 0 ? 0 : overlapIn;
    const retainEnd = N - nextOverlap;
    if (retainEnd > retainStart) {
      for (const f of clipFeatures(piece.doc.features, retainStart, retainEnd)) {
        doc.features.push({
          ...f,
          segments: f.segments.map((s) => ({ ...s, start: s.start + offset, end: s.end + offset })),
        });
      }
    }
    inputsMeta.push({
      index,
      name: piece.source?.name ?? piece.doc.name,
      seqLen: piece.source?.seq.length ?? N,
      circular: !!piece.source?.circular,
      role: i === 0 ? 'backbone' : 'insert',
      orientation: piece.reverse ? 'reverse' : 'forward',
      fragment: { start: retainStart, end: retainEnd },
      piece: { start: offset + retainStart, end: offset + retainEnd, len: Math.max(0, retainEnd - retainStart) },
      ends: { left: piece.left, right: piece.right },
      overlapIn: i === 0 ? 0 : overlapIn,
      overlapOut: nextOverlap,
      features: structuredClone(piece.source?.features ?? []),
    });
    if (i > 0) {
      const at = offset + overlapIn;
      const overlap = piece.doc.seq.slice(0, overlapIn);
      const ctx = junctionContext(productSeq, circular, at, FLANK);
      junctions.push({
        position: at,
        type: 'overlap',
        overhang: overlap,
        enzymes: [],
        leftInput: placements[i - 1].index,
        rightInput: index,
        sequence: { top: ctx.top, bottom: revcomp(ctx.top), start: ctx.start },
        regeneratedSites: [],
      });
    }
  }

  if (circular) {
    const overlap = best.acc.slice(best.acc.length - best.closing);
    const ctx = junctionContext(productSeq, circular, productLen, FLANK);
    junctions.push({
      position: productLen,
      type: 'overlap',
      overhang: overlap,
      enzymes: [],
      leftInput: placements[n - 1].index,
      rightInput: placements[0].index,
      sequence: { top: ctx.top, bottom: revcomp(ctx.top), start: ctx.start },
      regeneratedSites: [],
    });
  }

  doc.cloning = {
    type: 'gibson',
    operation: 'gibsonAssembly',
    created: new Date().toISOString(),
    circular,
    product: { name: doc.name, seqLen: productLen },
    inputs: inputsMeta,
    junctions,
  };

  const labels = inputsMeta.map((i) => `${i.name}${i.orientation === 'reverse' ? ' (reverse)' : ''}`);
  const overlaps = [...inputsMeta.slice(1).map((i) => i.overlapIn), ...(circular ? [best.closing] : [])];
  addHistory(
    doc,
    'Gibson assembly',
    `${labels.join(' + ')} → ${productLen} bp ${circular ? 'circular' : 'linear'}; overlaps ${overlaps.join(', ')} bp`,
  );
  doc.notes.description = description || provenanceText(inputsMeta, 'gibson', circular);

  const warnings = [];
  const cycles = new Set(finished.map((s) => cycleKey(s, circular)));
  if (cycles.size > 1) {
    warnings.push(`${cycles.size} possible assemblies; used the one with the largest total overlap.`);
  }
  if (inputsMeta.some((i) => i.orientation === 'reverse')) {
    warnings.push('At least one fragment was reverse-complemented to match the overlaps.');
  }
  return { doc, junctions, inputs: inputsMeta, order: placements.map((p) => ({ index: p.index, reverse: p.piece.reverse })), warnings };
}

// ---- history ---------------------------------------------------------------

function segmentRange(s) {
  return `${s.start + 1}-${s.end}`;
}

function richText(value) {
  const s = String(value ?? '');
  if (!s) return s;
  if (/^<html[\s>]/i.test(s)) return s;
  return `<html><body>${s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/\n/g, '<br>')}</body></html>`;
}

function featureNode(f) {
  const attrs = { ...f.extra, name: f.name, type: f.type };
  if (f.strand === 1) attrs.directionality = '1';
  else if (f.strand === -1) attrs.directionality = '2';
  const segs = (f.segments ?? []).map((s) =>
    el('Segment', { ...s.extra, range: segmentRange(s), color: s.color || f.color, name: s.name, type: s.extra?.type ?? 'standard' }),
  );
  const quals = (f.qualifiers ?? []).map((q) =>
    el(
      'Q',
      { name: q.name },
      (q.values ?? []).map((v) =>
        v.kind === 'int' ? el('V', { int: v.value }) : v.kind === 'predef' ? el('V', { predef: v.value }) : el('V', { text: richText(v.value) }),
      ),
    ),
  );
  return el('Feature', attrs, [...segs, ...quals]);
}

const featuresNode = (features) => el('Features', {}, (features ?? []).map(featureNode));

function rangeAttr(start, end) {
  return `${start}..${Math.max(start, end - 1)}`;
}

// Render doc.cloning as the binary format's history-tree XML (the content of its
// 0x07 segment, which the format stores compressed). Returns null when the
// document has no cloning record.
export function historyTreeXML(doc) {
  const c = doc?.cloning;
  if (!c) return null;
  const kids = [];

  const seenSites = new Set();
  for (const site of c.junctions.flatMap((j) => j.regeneratedSites ?? [])) {
    const key = `${site.name}:${site.pos}`;
    if (seenSites.has(key)) continue;
    seenSites.add(key);
    kids.push(el('RegeneratedSite', { name: site.name, pos: site.pos, siteCount: site.siteCount ?? 1 }));
  }

  // Display-only: the file format trims the overhang bases from these ranges;
  // the pieces' full extents are accurate enough for recolouring a map.
  const ranges = c.inputs.map((inp, i) =>
    el('ColorRange', { range: rangeAttr(inp.piece.start, inp.piece.end), colors: i === 0 ? 'black' : 'red' }),
  );
  kids.push(el('HistoryColors', {}, [el('StrandColors', { type: 'Product' }, [el('TopStrand', {}, ranges)])]));

  c.inputs.forEach((inp, i) => {
    let manipulation;
    if (c.type === 'gibson') manipulation = 'insert';
    else if (i === 0) manipulation = 'replace';
    else manipulation = inp.orientation === 'reverse' ? 'flipAndInsert' : 'insert';
    const left = inp.ends?.left;
    const right = inp.ends?.right;
    const hasEnzyme = left?.enzyme || right?.enzyme;
    kids.push(
      el('InputSummary', {
        manipulation,
        name1: left?.enzyme ?? undefined,
        name2: right?.enzyme ?? undefined,
        val1: inp.fragment ? inp.fragment.start : inp.piece.start,
        val2: inp.fragment ? inp.fragment.end : inp.piece.end,
        siteCount1: hasEnzyme ? 1 : undefined,
        siteCount2: hasEnzyme ? 1 : undefined,
      }),
    );
  });

  for (let i = 0; i < c.inputs.length; i++) {
    const inp = c.inputs[i];
    kids.push(
      el(
        'Node',
        {
          name: inp.name,
          type: 'DNA',
          seqLen: inp.seqLen,
          strandedness: 'double',
          ID: i + 1,
          circular: inp.circular ? '1' : undefined,
          operation: 'invalid',
          resurrectable: '1',
        },
        [featuresNode(inp.features)],
      ),
    );
  }

  const root = el(
    'Node',
    {
      name: doc.name,
      type: 'DNA',
      seqLen: doc.seq.length,
      strandedness: 'double',
      ID: 0,
      circular: doc.circular ? '1' : undefined,
      operation: c.operation,
    },
    kids,
  );
  return toXML(el('HistoryTree', {}, [root]), { declaration: true });
}

// Human-oriented summary for the command line or a panel.
export function cloningSummary(doc) {
  const c = doc?.cloning;
  if (!c) return null;
  return {
    type: c.type,
    operation: c.operation,
    size: doc.seq.length,
    circular: doc.circular,
    inputs: c.inputs.map((i) => ({
      name: i.name,
      role: i.role,
      size: i.seqLen,
      orientation: i.orientation,
      region: i.fragment,
    })),
    junctions: c.junctions.map((j) => ({
      position: j.position,
      kind: j.type === 'overlap' ? `${j.overhang.length} bp overlap` : j.type === 'blunt' ? 'blunt' : `${j.type}′ ${j.overhang}`,
      enzymes: (j.enzymes ?? []).filter(Boolean),
      regenerated: (j.regeneratedSites ?? []).map((s) => `${s.name} at ${s.pos + 1}`),
      sequence: j.sequence?.top ?? '',
    })),
  };
}
