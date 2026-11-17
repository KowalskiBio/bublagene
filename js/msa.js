// Multiple sequence alignment by the centre-star method: every sequence is
// aligned to one centre sequence (the longest), and the pairwise alignments are
// merged column by column so insertions relative to the centre are kept.
// Fast and good for similar sequences (variants, clones, reads of one locus);
// not a substitute for a progressive aligner on distant ones.

import { alignPair } from './align.js';

const GAP = '-';

// Merge a pairwise alignment (centre-with-gaps `pc`, other `ps`) into rows whose
// first row is the centre-with-gaps so far.
function merge(rows, pc, ps) {
  const out = rows.map(() => '');
  let added = '';
  const A = rows[0].length;
  const B = pc.length;
  let a = 0;
  let b = 0;
  while (a < A || b < B) {
    const ca = a < A ? rows[0][a] : null;
    const cb = b < B ? pc[b] : null;
    if (cb === null || (ca === GAP && cb !== GAP)) {          // a column only the earlier alignment has
      rows.forEach((r, i) => { out[i] += r[a]; });
      added += GAP;
      a++;
    } else if (ca === null || (cb === GAP && ca !== GAP)) {   // a new insertion
      rows.forEach((_, i) => { out[i] += GAP; });
      added += ps[b];
      b++;
    } else {                                                  // the same column in both
      rows.forEach((r, i) => { out[i] += r[a]; });
      added += ps[b];
      a++;
      b++;
    }
  }
  return [...out, added];
}

/**
 * @param {{name:string, seq:string}[]} items
 * @returns {{names:string[], rows:string[], consensus:string, conservation:number[], centre:number, identity:number}}
 */
export function multipleAlign(items, { maxLength = 8000, ...alignOptions } = {}) {
  if (items.length < 2) throw new Error('Pick at least two sequences to align');
  for (const it of items) {
    if (it.seq.length > maxLength) throw new Error(`${it.name} is ${it.seq.length} bp; multiple alignment is limited to ${maxLength} bp each`);
  }
  const seqs = items.map((it) => it.seq.toUpperCase());
  const centre = seqs.reduce((best, s, i) => (s.length > seqs[best].length ? i : best), 0);
  let rows = [seqs[centre]];
  const order = [centre];
  for (let i = 0; i < seqs.length; i++) {
    if (i === centre) continue;
    const pair = alignPair(seqs[centre], seqs[i], alignOptions);
    rows = merge(rows, pair.alignedA, pair.alignedB);
    order.push(i);
  }
  // put the rows back in the order the sequences were given
  const byInput = new Array(items.length);
  order.forEach((inputIndex, k) => { byInput[inputIndex] = rows[k]; });
  const width = byInput[0].length;
  let consensus = '';
  const conservation = [];
  let same = 0;
  for (let c = 0; c < width; c++) {
    const counts = {};
    for (const r of byInput) counts[r[c]] = (counts[r[c]] ?? 0) + 1;
    const [top, n] = Object.entries(counts).sort((x, y) => y[1] - x[1])[0];
    conservation.push(n / byInput.length);
    consensus += n / byInput.length > 0.5 ? top : GAP;
    if (n === byInput.length && top !== GAP) same++;
  }
  return { names: items.map((it) => it.name), rows: byInput, consensus, conservation, centre, identity: width ? (100 * same) / width : 0 };
}

// Columns where the rows differ, with positions counted along the first row.
export function alignmentVariants(rows) {
  const out = [];
  let refPos = 0;
  for (let c = 0; c < rows[0].length; c++) {
    const ref = rows[0][c];
    const bases = rows.map((r) => r[c]);
    if (new Set(bases).size > 1) {
      out.push({ column: c, position: refPos, ref, bases, type: ref === GAP ? 'insertion' : bases.includes(GAP) ? 'deletion' : 'substitution' });
    }
    if (ref !== GAP) refPos++;
  }
  return out;
}
