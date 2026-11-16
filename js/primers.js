// Primer binding, design and PCR simulation.

import { revcomp, meltingTemp, gcContent, mod, subseq } from './seq.js';
import { newDocument, newId, addHistory, featuresInRegion } from './model.js';

// Where does a primer anneal? Requires an exact match of the 3'-terminal
// `minAnneal` bases, then extends towards the 5' end while bases still match.
// Returns [{strand, start, end, annealed, tm, tail}] in top-strand coordinates.
export function findBindingSites(seq, circular, primer, minAnneal = 15) {
  const p = primer.toUpperCase();
  const L = seq.length;
  if (p.length < minAnneal || L < minAnneal) return [];
  const hay = circular ? seq + seq.slice(0, p.length) : seq;
  const sites = [];
  // forward: primer 3' end = top strand, read left to right
  const core = p.slice(-minAnneal);
  for (let i = hay.indexOf(core); i >= 0 && i < L + (circular ? 0 : 0); i = hay.indexOf(core, i + 1)) {
    const end3 = i + minAnneal; // exclusive, top coords (may exceed L when circular)
    let a = minAnneal;
    while (a < p.length) {
      const pos = end3 - a - 1;
      const base = circular ? seq[mod(pos, L)] : seq[pos];
      if (pos < 0 && !circular) break;
      if (base !== p[p.length - a - 1]) break;
      a++;
    }
    const start = end3 - a;
    sites.push(mk(1, start, end3, p.slice(p.length - a), p.slice(0, p.length - a)));
    if (i >= L - 1) break;
  }
  // reverse: primer anneals to the top strand as its reverse complement
  const rc = revcomp(p);
  const coreR = rc.slice(0, minAnneal);
  for (let i = hay.indexOf(coreR); i >= 0 && i < L; i = hay.indexOf(coreR, i + 1)) {
    let a = minAnneal;
    while (a < p.length) {
      const pos = i + a;
      if (!circular && pos >= L) break;
      if (seq[mod(pos, L)] !== rc[a]) break;
      a++;
    }
    sites.push(mk(-1, i, i + a, p.slice(p.length - a), p.slice(0, p.length - a)));
  }
  function mk(strand, start, end, annealed, tail) {
    return {
      strand,
      start: circular ? mod(start, L) : start,
      end: circular ? mod(end, L) || L : end,
      annealed,
      tail,
      tm: meltingTemp(annealed),
    };
  }
  // de-duplicate (the same site can be found twice through the wrap-around copy)
  const seen = new Set();
  return sites.filter((s) => {
    const k = `${s.strand}:${s.start}`;
    if (seen.has(k)) return false;
    seen.add(k);
    return true;
  });
}

// Pick the primer of length minLen..maxLen whose Tm is closest to the target.
function pickLength(template, targetTm, minLen, maxLen) {
  let best = null;
  for (let n = minLen; n <= Math.min(maxLen, template.length); n++) {
    const s = template.slice(0, n);
    const tm = meltingTemp(s);
    // small bonus for a G/C clamp at the 3' end
    const clamp = /[GC]$/.test(s) ? 0 : 1.5;
    const score = Math.abs(tm - targetTm) + clamp;
    if (!best || score < best.score) best = { seq: s, tm, score };
  }
  return best;
}

// Design a forward/reverse primer pair that amplifies [start, end).
export function designPrimers(seq, circular, start, end, { targetTm = 60, minLen = 18, maxLen = 32, fwdTail = '', revTail = '' } = {}) {
  const region = subseq(seq, start, end, circular);
  if (region.length < 2 * minLen) throw new Error(`Region too short (${region.length} bp) for two ${minLen}-nt primers`);
  const f = pickLength(region, targetTm, minLen, maxLen);
  const r = pickLength(revcomp(region), targetTm, minLen, maxLen);
  const mkp = (core, tail, dir) => ({
    seq: (tail + core.seq).toUpperCase(),
    annealing: core.seq,
    tm: core.tm,
    gc: gcContent(core.seq),
    dir,
  });
  return { forward: mkp(f, fwdTail, 1), reverse: mkp(r, revTail, -1), productLength: region.length + fwdTail.length + revTail.length };
}

// Simulate PCR on a template document with two primers (sequences 5'->3').
export function simulatePCR(doc, fwdSeq, revSeq, { minAnneal = 15 } = {}) {
  const L = doc.seq.length;
  const fSites = findBindingSites(doc.seq, doc.circular, fwdSeq, minAnneal);
  const rSites = findBindingSites(doc.seq, doc.circular, revSeq, minAnneal);
  // allow either primer to act as forward
  const pairs = [];
  for (const [aSeq, bSeq, aS, bS] of [[fwdSeq, revSeq, fSites, rSites], [revSeq, fwdSeq, rSites, fSites]]) {
    for (const a of aS.filter((s) => s.strand === 1)) {
      for (const b of bS.filter((s) => s.strand === -1)) {
        const inner = doc.circular ? mod(b.end - a.start, L) || L : b.end - a.start;
        if (inner <= 0 || (!doc.circular && b.end <= a.start)) continue;
        pairs.push({ a, b, aSeq, bSeq, inner });
      }
    }
  }
  if (!pairs.length) throw new Error('No product: the primers do not bind facing each other on this template (need ≥15 matching bases at each 3′ end).');
  pairs.sort((x, y) => x.inner - y.inner);
  const { a, b, aSeq, bSeq } = pairs[0];
  const templatePart = subseq(doc.seq, a.start, b.end, doc.circular);
  const product = (a.tail + templatePart + revcomp(b.tail)).toUpperCase();
  const out = newDocument({ name: `${doc.name} PCR product`, seq: product, circular: false });
  out.features = featuresInRegion(doc, a.start, templatePart.length, a.tail.length);
  out.primers = [
    { id: newId(), name: 'forward primer', seq: aSeq.toUpperCase(), description: '', extra: {} },
    { id: newId(), name: 'reverse primer', seq: bSeq.toUpperCase(), description: '', extra: {} },
  ];
  out.history = [...(doc.history || [])];
  addHistory(out, 'PCR', `PCR from ${doc.name} (${a.start + 1}–${b.end}), ${product.length} bp`);
  out.notes.description = `PCR product from ${doc.name}`;
  return { doc: out, warnings: pairs.length > 1 ? [`${pairs.length} possible products; the shortest (${product.length} bp) was used.`] : [] };
}
