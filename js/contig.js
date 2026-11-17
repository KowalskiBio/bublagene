// Contig assembly from Sanger-style reads.
//
// Overlaps are suffix/prefix matches scored by length and mismatch rate.
// Bases are merged by quality-weighted voting; a tie between different bases
// becomes N. Reads may be given in either orientation.

import { revcomp } from './seq.js';

const DEFAULT_Q = 20;

function normaliseRead(read, index) {
  if (typeof read === 'string') return { name: `read ${index + 1}`, seq: read.toUpperCase(), qualities: null };
  const seq = String(read.seq ?? '').toUpperCase();
  return { name: read.name ?? `read ${index + 1}`, seq, qualities: read.qualities ?? read.quality ?? null };
}

function qAt(read, i) {
  const q = read.qualities;
  if (!q) return DEFAULT_Q;
  if (Array.isArray(q)) return Number(q[i] ?? DEFAULT_Q);
  return DEFAULT_Q;
}

function complementBase(b) {
  return { A: 'T', C: 'G', G: 'C', T: 'A', N: 'N' }[b] ?? 'N';
}

// Best suffix(A) / prefix(B) overlap. `reverseB` is handled by the caller.
export function overlapScore(a, b, { minOverlap = 30, maxMismatchRate = 0.1 } = {}) {
  const max = Math.min(a.length, b.length);
  let best = null;
  for (let k = max; k >= minOverlap; k--) {
    let mismatches = 0;
    for (let i = 0; i < k; i++) {
      const x = a[a.length - k + i];
      const y = b[i];
      if (x !== y) {
        mismatches++;
      }
    }
    const allowed = Math.floor(maxMismatchRate * k);
    if (mismatches <= allowed) {
      const score = k - 2 * mismatches;
      if (!best || score > best.score) best = { length: k, mismatches, score };
    }
  }
  return best;
}

function columnsFor(read, reverse = false) {
  const seq = reverse ? revcomp(read.seq) : read.seq;
  return [...seq].map((base, i) => {
    const q = qAt(read, reverse ? read.seq.length - 1 - i : i);
    return { votes: { [base]: q }, depth: 1, base, q };
  });
}

function emptyColumns(length) {
  return Array.from({ length }, () => ({ votes: {}, depth: 0, base: 'N', q: 0 }));
}

function callColumn(votes, depth) {
  let best = null;
  let bestQ = -1;
  let tie = false;
  for (const [base, q] of Object.entries(votes)) {
    if (q > bestQ) { best = base; bestQ = q; tie = false; }
    else if (q === bestQ && base !== best) tie = true;
  }
  return { votes, depth, base: tie ? 'N' : best ?? 'N', q: tie ? 0 : bestQ };
}

function mergeColumns(left, right) {
  const out = emptyColumns(Math.max(left.length, right.length));
  for (let i = 0; i < out.length; i++) {
    const votes = {};
    let depth = 0;
    for (const side of [left, right]) {
      const col = side[i];
      if (!col) continue;
      depth += col.depth;
      for (const [base, q] of Object.entries(col.votes)) votes[base] = (votes[base] ?? 0) + q;
    }
    out[i] = callColumn(votes, depth);
  }
  return out;
}

function mergeContigs(A, B, ov, reverseB) {
  const bSeq = reverseB ? revcomp(B.seq) : B.seq;
  const bColumns = reverseB ? B.columns.map((c) => {
    const votes = {};
    for (const [base, q] of Object.entries(c.votes)) votes[complementBase(base)] = q;
    return callColumn(votes, c.depth);
  }).reverse() : B.columns;
  const offset = A.seq.length - ov;
  const length = Math.max(A.seq.length, offset + bSeq.length);
  const left = [...A.columns, ...emptyColumns(Math.max(0, length - A.columns.length))];
  const right = emptyColumns(length);
  for (let i = 0; i < bColumns.length; i++) right[offset + i] = bColumns[i];
  const columns = mergeColumns(left, right);
  const seq = columns.map((c) => c.base).join('');
  const layout = [
    ...A.layout,
    ...B.layout.map((l) => ({
      name: l.name,
      offset: offset + (reverseB ? bSeq.length - l.offset - l.length : l.offset),
      length: l.length,
      reverse: reverseB ? !l.reverse : l.reverse,
    })),
  ].sort((x, y) => x.offset - y.offset);
  return {
    seq,
    columns,
    layout,
    reads: [...A.reads, ...B.reads],
    length: seq.length,
    depth: columns.reduce((n, c) => n + c.depth, 0) / Math.max(1, columns.length),
  };
}

export function assembleContigs(reads, {
  minOverlap = 30, maxMismatchRate = 0.1, maxContigs = 1000,
} = {}) {
  const list = reads.map(normaliseRead).filter((r) => r.seq.length);
  if (!list.length) return [];
  let contigs = list.map((read) => ({
    seq: read.seq,
    columns: columnsFor(read),
    layout: [{ name: read.name, offset: 0, length: read.seq.length, reverse: false }],
    reads: [read.name],
    length: read.seq.length,
    depth: 1,
  }));

  let changed = true;
  while (changed && contigs.length > 1) {
    changed = false;
    let best = null;
    for (let i = 0; i < contigs.length; i++) {
      for (let j = 0; j < contigs.length; j++) {
        if (i === j) continue;
        for (const reverse of [false, true]) {
          const B = reverse ? { ...contigs[j], seq: revcomp(contigs[j].seq) } : contigs[j];
          const ov = overlapScore(contigs[i].seq, B.seq, { minOverlap, maxMismatchRate });
          if (ov && (!best || ov.score > best.ov.score)) best = { i, j, reverse, ov };
        }
      }
    }
    if (best) {
      const merged = mergeContigs(contigs[best.i], contigs[best.j], best.ov.length, best.reverse);
      const next = contigs.filter((_, k) => k !== best.i && k !== best.j);
      next.push(merged);
      contigs = next;
      changed = true;
      if (contigs.length > maxContigs) break;
    }
  }
  return contigs
    .map((c) => ({
      seq: c.seq,
      length: c.length,
      reads: c.reads,
      layout: c.layout,
      depth: Math.round(c.depth * 100) / 100,
      circular: false,
    }))
    .sort((a, b) => b.length - a.length);
}

// Majority/quality consensus of reads already placed at known offsets.
export function consensusAt(reference, reads, offset = 0) {
  const columns = emptyColumns(reference.length);
  for (const read of reads.map(normaliseRead)) {
    for (let i = 0; i < read.seq.length; i++) {
      const p = i + offset;
      if (p < 0 || p >= columns.length) continue;
      const base = read.seq[i];
      const q = qAt(read, i);
      columns[p].votes[base] = (columns[p].votes[base] ?? 0) + q;
      columns[p].depth++;
    }
  }
  return columns.map((c) => callColumn(c.votes, c.depth).base).join('');
}
