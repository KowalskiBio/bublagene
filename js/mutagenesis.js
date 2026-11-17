// Site-directed mutagenesis: design the primers for a substitution, insertion
// or deletion in a plasmid, and build the expected mutant.
//
// Two designs are offered, both amplifying the whole circular template:
//  - 'overlap'       a complementary pair centred on the change (QuikChange style).
//                    Flanks grow until the classic estimate
//                    Tm = 81.5 + 0.41·GC% − 675/N − 100·mismatch/N reaches the target.
//  - 'backtoback'    two primers that point away from each other, meeting at the
//                    change (Q5-SDM style). Only the annealing parts are sized
//                    to the target Tm; the amplicon is then phosphorylated and ligated.

import { revcomp, meltingTemp, gcContent, mod, rangeLength } from './seq.js';
import { cloneDoc, replaceRange, addHistory, newId, frameOffset, isCodingSegment, qualifierText } from './model.js';
import { getCodonTable } from './codon-tables.js';

const around = (seq, from, len) => {
  const L = seq.length;
  const start = mod(from, L);
  return (seq + seq).slice(start, start + len);
};

// ---- the change --------------------------------------------------------------

// Where residue `residue` (1-based) of a CDS lives in the sequence, and what to
// write there to make it `newAa`. The codon needing the fewest base changes wins.
export function aminoAcidChange(doc, feature, residue, newAa) {
  if (feature.type !== 'CDS') throw new Error('Pick a CDS feature');
  const L = doc.seq.length;
  const order = [];
  for (const seg of feature.segments.filter(isCodingSegment)) {
    for (let k = 0; k < rangeLength(seg, L); k++) order.push((seg.start + k) % L);
  }
  if (feature.strand === -1) order.reverse();
  const off = frameOffset(feature);
  const at = off + (residue - 1) * 3;
  if (!(residue >= 1) || at + 3 > order.length) throw new Error(`${feature.name} has no residue ${residue}`);
  const pos = [order[at], order[at + 1], order[at + 2]];
  const step = feature.strand === -1 ? -1 : 1;
  if (pos[1] !== mod(pos[0] + step, L) || pos[2] !== mod(pos[1] + step, L)) throw new Error('That codon is split across segments; edit its bases directly');
  const table = getCodonTable(Number(qualifierText(feature, 'transl_table')) || 1);
  const start = feature.strand === -1 ? pos[2] : pos[0];
  const topOld = around(doc.seq, start, 3);
  const oldCodon = feature.strand === -1 ? revcomp(topOld) : topOld;
  const oldAa = table.codons[oldCodon] ?? 'X';
  const aa = String(newAa).toUpperCase() === 'STOP' ? '*' : String(newAa).toUpperCase();
  const options = Object.keys(table.codons).filter((c) => table.codons[c] === aa);
  if (!options.length) throw new Error(`No codon for “${newAa}”`);
  const dist = (c) => [...c].filter((ch, i) => ch !== oldCodon[i]).length;
  const newCodon = [...options].sort((a, b) => dist(a) - dist(b) || a.localeCompare(b))[0];
  if (newCodon === oldCodon) throw new Error(`Residue ${residue} is already ${aa === '*' ? 'a stop' : aa}`);
  return {
    start, end: start + 3 > L ? start + 3 - L : start + 3,
    oldBases: topOld, newBases: feature.strand === -1 ? revcomp(newCodon) : newCodon,
    oldCodon, newCodon, oldAa, newAa: aa,
    label: `${oldAa}${residue}${aa === '*' ? '*' : aa}`,
  };
}

// ---- primers -----------------------------------------------------------------

const classicTm = (seq, mismatches) => {
  const n = seq.length;
  return 81.5 + 0.41 * gcContent(seq) - 675 / n - (100 * mismatches) / n;
};

function overlapPair(template, start, end, mut, { targetTm, maxLen }) {
  const L = template.length;
  const delLen = end - start;
  const mismatches = Math.max(mut.length, delLen);
  for (let flank = 10; flank <= 40; flank++) {
    const left = around(template, start - flank, flank);
    const right = around(template, end, flank);
    const seq = left + mut + right;
    if (seq.length > maxLen) break;
    if (classicTm(seq, mismatches) >= targetTm) return { seq, left: flank, right: flank, tm: classicTm(seq, mismatches) };
  }
  // could not reach the target within the length cap: take the longest allowed
  const flank = Math.max(10, Math.floor((maxLen - mut.length) / 2));
  const seq = around(template, start - flank, flank) + mut + around(template, end, flank);
  return { seq, left: flank, right: flank, tm: classicTm(seq, mismatches), short: true };
}

function annealing(template, from, dir, targetTm, { min = 15, max = 40 } = {}) {
  // grow from `from` outward (dir +1 reads right, -1 reads left) until the Tm is met
  for (let n = min; n <= max; n++) {
    const part = dir > 0 ? around(template, from, n) : around(template, from - n, n);
    if (meltingTemp(part) >= targetTm) return { part, n, tm: meltingTemp(part) };
  }
  const part = dir > 0 ? around(template, from, max) : around(template, from - max, max);
  return { part, n: max, tm: meltingTemp(part), short: true };
}

/**
 * @param {any} doc a circular plasmid
 * @param {{start:number,end:number,newBases:string,label?:string}} change bases [start,end) are replaced by newBases
 */
export function designMutagenesis(doc, change, { method = 'overlap', targetTm = null, maxLen = 48 } = {}) {
  if (!doc.circular) throw new Error('Mutagenesis primers amplify the whole plasmid, so the template must be circular');
  const L = doc.seq.length;
  const { start, end } = change;
  const mut = String(change.newBases ?? '').toUpperCase().replace(/[^ACGT]/g, '');
  if (!(start >= 0 && end >= start && end <= L)) throw new Error('The change must sit inside the sequence without wrapping the origin');
  if (!mut && end === start) throw new Error('Nothing to change');
  const old = doc.seq.slice(start, end);
  if (mut === old) throw new Error('The new bases are the same as the old ones');
  const kind = !old ? 'insertion' : !mut ? 'deletion' : old.length === mut.length ? 'substitution' : 'replacement';
  const notes = [];
  let forward;
  let reverse;
  let tmTarget = targetTm ?? (method === 'overlap' ? 78 : 62);

  if (method === 'overlap') {
    if (mut.length > maxLen - 20) throw new Error(`An insertion of ${mut.length} bp is too long for one overlapping pair; use back-to-back primers`);
    const p = overlapPair(doc.seq, start, end, mut, { targetTm: tmTarget, maxLen });
    forward = { seq: p.seq, tm: p.tm };
    reverse = { seq: revcomp(p.seq), tm: p.tm };
    if (p.short) notes.push(`Could not reach Tm ${tmTarget} °C within ${maxLen} nt; the primers are as long as allowed.`);
    notes.push(`Primers overlap completely and carry the change in the middle, with ${p.left} nt of template on each side.`);
    notes.push('Amplify the whole plasmid, digest the template with DpnI, then transform.');
  } else if (method === 'backtoback') {
    const half = Math.ceil(mut.length / 2);
    const long = mut.length > 18;
    const fTail = long ? mut.slice(half) : mut;
    const rTail = long ? revcomp(mut.slice(0, half)) : '';
    const f = annealing(doc.seq, end, 1, tmTarget);
    const r = annealing(doc.seq, start, -1, tmTarget);
    forward = { seq: fTail + f.part, tm: f.tm, tail: fTail, anneal: f.part };
    reverse = { seq: rTail + revcomp(r.part), tm: r.tm, tail: rTail, anneal: revcomp(r.part) };
    if (f.short || r.short) notes.push(`Could not reach Tm ${tmTarget} °C with a 40 nt annealing part.`);
    notes.push('Primers point away from each other and meet at the change; the amplicon is the whole plasmid.');
    notes.push('Phosphorylate and ligate the product (KLD or T4 PNK + T4 ligase), digest the template, then transform.');
  } else {
    throw new Error(`Unknown method ${method}`);
  }

  const mutant = cloneDoc(doc);
  mutant.id = newId();
  replaceRange(mutant, start, end, mut);
  const label = change.label ?? `${kind} at ${start + 1}`;
  mutant.name = `${doc.name} ${label}`;
  mutant.source = { format: 'new' };
  const name = (dir) => `${label} ${dir}`;
  forward = { ...forward, name: name('fwd'), length: forward.seq.length, gc: gcContent(forward.seq) };
  reverse = { ...reverse, name: name('rev'), length: reverse.seq.length, gc: gcContent(reverse.seq) };
  mutant.primers = [
    ...mutant.primers,
    { id: newId(), name: forward.name, seq: forward.seq, description: `Mutagenic primer for ${label}`, extra: {} },
    { id: newId(), name: reverse.name, seq: reverse.seq, description: `Mutagenic primer for ${label}`, extra: {} },
  ];
  addHistory(mutant, 'Mutagenesis', `${label}: ${old || '–'} → ${mut || '–'} at ${start + 1}–${end || start} (${method})`);
  return { kind, method, label, targetTm: tmTarget, forward, reverse, mutant, notes, old, mut, delta: mut.length - old.length };
}
