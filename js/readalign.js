// Place a Sanger read on a reference.
//
// The read is seeded against the reference with shared k-mers in both
// orientations, which finds the right stretch of a long sequence cheaply. Only
// that stretch (plus some slack for indels) is aligned in full, so a 1 kb read
// against a 100 kb plasmid costs about the same as against a 2 kb one.

import { alignLocal } from './align.js';
import { revcomp } from './seq.js';

const K = 10;
const BIN = 32;      // diagonals this close count as the same placement
const SLACK = 80;    // reference either side of the seeded stretch

// Reference k-mer index: k-mer -> start positions.
function indexKmers(seq) {
  const index = new Map();
  for (let i = 0; i + K <= seq.length; i++) {
    const kmer = seq.slice(i, i + K);
    if (/[^ACGT]/.test(kmer)) continue;
    const list = index.get(kmer);
    if (list) list.push(i); else index.set(kmer, [i]);
  }
  return index;
}

// The best placement of `read` on the indexed reference: how many k-mers agree
// on one diagonal band, and the middle of that band.
function seed(index, read) {
  const votes = new Map();
  for (let q = 0; q + K <= read.length; q++) {
    const hits = index.get(read.slice(q, q + K));
    if (!hits || hits.length > 20) continue; // repeats say little about where the read sits
    for (const p of hits) {
      const bin = Math.floor((p - q) / BIN);
      votes.set(bin, (votes.get(bin) ?? 0) + 1);
    }
  }
  let best = null;
  for (const [bin, n] of votes) {
    // a band is a pair of neighbouring bins, so a placement on a boundary is not split
    const score = n + (votes.get(bin + 1) ?? 0);
    if (!best || score > best.score) best = { score, bin };
  }
  return best ? { score: best.score, diagonal: (best.bin + 1) * BIN } : null;
}

/**
 * Turn an aligned pair of strings into a row the sequence and alignment views
 * draw: one character per reference position the read covers ('-' where the
 * read lacks a base), insertions keyed by the position they sit before, and
 * `map`, the reference position of each read base (null outside the alignment).
 * `n` is the length of the read, in the orientation that was aligned.
 */
export function rowFromAlignment({ alignedRef, alignedRead, startRef, startRead }, n) {
  const chars = [];
  const insertions = new Map();
  const map = new Array(n).fill(null);
  let q = startRead;
  for (let i = 0; i < alignedRef.length; i++) {
    const a = alignedRef[i];
    const b = alignedRead[i];
    if (a !== '-' && b !== '-') { map[q++] = startRef + chars.length; chars.push(b); }
    else if (a !== '-') chars.push('-');
    else {
      const at = startRef + chars.length;
      insertions.set(at, (insertions.get(at) ?? '') + b);
      map[q++] = at - 1 >= startRef ? at - 1 : null;
    }
  }
  return { start: startRef, chars, insertions, map, end: startRef + chars.length };
}

/**
 * Align a trace to a reference.
 * Returns null when it does not match. Otherwise a read in the shape the
 * sequence view draws (start, chars, insertions, label, reverse), plus `map`
 * (the reference position of each called base, or null) and the identity.
 */
export function placeRead(reference, trace, { label = trace.name || 'read', minIdentity = 70, minAligned = 30 } = {}) {
  const ref = String(reference).toUpperCase();
  const bases = String(trace.bases ?? '').toUpperCase();
  if (ref.length < K || bases.length < K) return null;
  const index = indexKmers(ref);
  const fwd = seed(index, bases);
  const rev = seed(index, revcomp(bases));
  const reverse = !!rev && (!fwd || rev.score > fwd.score);
  const hit = reverse ? rev : fwd;
  if (!hit || hit.score < 4) return null;
  const oriented = reverse ? revcomp(bases) : bases;

  const from = Math.max(0, hit.diagonal - SLACK);
  const to = Math.min(ref.length, hit.diagonal + oriented.length + SLACK + BIN);
  if (to - from < K) return null;
  if ((to - from) * oriented.length > 8_000_000) return null; // too large to align in full
  const aln = alignLocal(ref.slice(from, to), oriented);
  const aligned = aln.matches + aln.mismatches;
  if (aligned < minAligned || aln.identityAligned < minIdentity) return null;

  const n = oriented.length;
  const row = rowFromAlignment(
    { alignedRef: aln.alignedA, alignedRead: aln.alignedB, startRef: from + aln.startA, startRead: aln.startB },
    n,
  );
  const { start, chars, insertions } = row;
  const map = reverse ? row.map.reverse() : row.map;
  return {
    label, reverse, start, chars, insertions, map,
    end: start + chars.length,
    identity: aln.identityAligned,
    matches: aln.matches, mismatches: aln.mismatches,
    readStart: reverse ? n - aln.endB : aln.startB,
    readEnd: reverse ? n - aln.startB : aln.endB,
  };
}

/**
 * The row for a pairwise alignment where the first sequence is the query and
 * the second is the reference (the Alignment panel's convention).
 */
export function rowFromPairwise(result, queryLength, label) {
  const row = rowFromAlignment(
    { alignedRef: result.alignedB, alignedRead: result.alignedA, startRef: result.startB, startRead: result.startA },
    queryLength,
  );
  return { label, reverse: false, ...row, identity: result.identityAligned, matches: result.matches, mismatches: result.mismatches };
}
