// Primer binding, design, thermodynamics and PCR simulation.
//
// Binding is mismatch-aware: the 3' end must match exactly for `threePrimeExact`
// bases (default 1) and at most `maxMismatches` differences are allowed in the
// annealed part. A 5' tail is carried into the product but never anneals.
//
// Design ranks candidate pairs by Tm match, Tm difference, GC, 3' clamp,
// secondary structure and how many extra binding sites each primer has.

import { revcomp, meltingTemp, gcContent, mod, subseq } from './seq.js';
import { primerProperties, duplexThermo, hairpin, dimer, pairDimer, clampScore } from './thermo.js';
import { newDocument, newId, addHistory, featuresInRegion } from './model.js';
import { ENZYME_BY_NAME, LADDERS, migration } from './enzymes.js';

const mismatch = (a, b) => (a === b ? 0 : 1);

// Where does a primer anneal? Returns [{strand, start, end, annealed, tail, tm,
// mismatches, exact}] in top-strand coordinates. `end` for a forward site is the
// 3' end; for a reverse site `start` is the 3' end (the primer reads right-to-left).
export function findBindingSites(seq, circular, primer, minAnneal = 15, {
  maxMismatches = 0, threePrimeExact = 1, maxSites = 200,
} = {}) {
  const p = String(primer).toUpperCase().replace(/U/g, 'T');
  const L = seq.length;
  if (p.length < minAnneal || L < minAnneal) return [];
  const hay = circular ? seq + seq.slice(0, Math.max(0, p.length - 1)) : seq;
  const sites = [];

  // Forward: the query's 3' end is on the right; a 5' tail is left of the
  // annealed part. Reverse: the query is revcomp(primer), so the primer's 3'
  // end is on the left and its 5' tail is right of the annealed part.
  const scan = (query, strand) => {
    const q = query;
    const n = q.length;
    const anchor = Math.min(threePrimeExact, n);
    const limitEnd = circular ? L + n : Math.min(hay.length, L + n);
    for (let pos = 0; pos + minAnneal <= limitEnd; pos++) {
      if (strand === 1) {
        const end = pos + minAnneal;
        if (end > hay.length) break;
        let exact = true;
        for (let k = 0; k < anchor; k++) if (hay[end - 1 - k] !== q[n - 1 - k]) { exact = false; break; }
        if (!exact) continue;
        let best = null;
        for (let a = Math.min(n, end); a >= minAnneal; a--) {
          let mm = 0;
          for (let k = 0; k < a; k++) if (hay[end - a + k] !== q[n - a + k]) mm++;
          if (mm <= maxMismatches) { best = { a, start: end - a, mm }; break; }
        }
        if (!best) continue;
        const annealed = q.slice(n - best.a);
        const tail = q.slice(0, n - best.a);
        push(strand, best.start, best.start + best.a, annealed, tail, best.mm);
      } else {
        const start = pos;
        if (start + minAnneal > hay.length) break;
        let exact = true;
        for (let k = 0; k < anchor; k++) if (hay[start + k] !== q[k]) { exact = false; break; }
        if (!exact) continue;
        let best = null;
        for (let a = Math.min(n, hay.length - start); a >= minAnneal; a--) {
          let mm = 0;
          for (let k = 0; k < a; k++) if (hay[start + k] !== q[k]) mm++;
          if (mm <= maxMismatches) { best = { a, start, mm }; break; }
        }
        if (!best) continue;
        const annealed = revcomp(q.slice(0, best.a));
        const tail = revcomp(q.slice(best.a));
        push(-1, best.start, best.start + best.a, annealed, tail, best.mm);
      }
    }
  };

  function push(strand, start, end, annealed, tail, mm) {
    if (sites.length >= maxSites) return;
    if (circular) {
      sites.push({
        strand, start: mod(start, L), end: mod(end, L) || L,
        annealed, tail, tm: meltingTemp(annealed), mismatches: mm,
        exact: mm === 0 && tail.length === 0,
      });
    } else {
      if (start < 0 || end > L) return;
      sites.push({
        strand, start, end, annealed, tail, tm: meltingTemp(annealed), mismatches: mm,
        exact: mm === 0 && tail.length === 0,
      });
    }
  }

  scan(p, 1);
  scan(revcomp(p), -1);

  const seen = new Set();
  return sites.filter((s) => {
    const k = `${s.strand}:${s.start}:${s.end}`;
    if (seen.has(k)) return false;
    seen.add(k);
    return true;
  });
}

// ---- design -----------------------------------------------------------------

function tailFor(spec, enzyme, extraClamp = 'GC') {
  const site = ENZYME_BY_NAME[spec];
  const seq = String(site ? site.site : spec ?? '');
  const clamp = seq ? extraClamp.slice(0, Math.max(0, 6 - seq.length)) : '';
  return (clamp + seq).toUpperCase();
}

function pickCandidate(template, targetTm, minLen, maxLen, options) {
  let best = null;
  for (let n = minLen; n <= Math.min(maxLen, template.length); n++) {
    const core = template.slice(0, n);
    if (/[^ACGT]/.test(core)) continue;
    const props = primerProperties(core, options);
    if (props.tm < (options.minTm ?? -Infinity) || props.tm > (options.maxTm ?? Infinity)) continue;
    const clampPenalty = props.clamp.threePrimeGc ? 0 : 1.5;
    const gcPenalty = props.gc < 35 || props.gc > 65 ? 4 : 0;
    const hairpinPenalty = props.hairpin.dg < -3 ? (-3 - props.hairpin.dg) * 2 : 0;
    const score = Math.abs(props.tm - targetTm) + clampPenalty + gcPenalty + hairpinPenalty;
    if (!best || score < best.score) best = { seq: core, tm: props.tm, props, score, penalty: clampPenalty + gcPenalty + hairpinPenalty };
  }
  return best;
}

// Design a forward/reverse pair for [start, end). Backwards-compatible wrapper.
export function designPrimers(seq, circular, start, end, opts = {}) {
  const result = designPrimerPairs(seq, circular, start, end, opts);
  return result;
}

export function designPrimerPairs(seq, circular, start, end, {
  targetTm = 60, minTm = 55, maxTm = 68, minLen = 18, maxLen = 32,
  maxTmDiff = 4, minProduct = 0, maxProduct = Infinity, flex = 0,
  fwdTail = '', revTail = '', fwdRestriction = null, revRestriction = null,
  maxCandidates = 24, specificity = true, minAnneal = 15, salt = {},
} = {}) {
  const L = seq.length;
  let start0 = mod(start, L);
  let end0 = end;
  if (!circular && end0 <= start0) end0 = L;
  const window = circular ? subseq(seq, start0, end0, true) : seq.slice(start0, end0);
  if (window.length < 2 * minLen) throw new Error(`Region too short (${window.length} bp) for two ${minLen}-nt primers`);

  const forwardTail = fwdTail || (fwdRestriction ? tailFor(fwdRestriction) : '');
  const reverseTail = revTail || (revRestriction ? tailFor(revRestriction) : '');
  const checks = [];
  const candidates = [];
  const shifts = [];
  for (let a = 0; a <= Math.min(flex, Math.max(0, window.length - 2 * minLen)); a++) {
    for (let b = 0; b <= Math.min(flex, Math.max(0, window.length - a - 2 * minLen)); b++) shifts.push([a, b]);
  }
  if (!shifts.length) shifts.push([0, 0]);

  for (const [leftTrim, rightTrim] of shifts) {
    const region = window.slice(leftTrim, window.length - rightTrim);
    const productLength = region.length + forwardTail.length + reverseTail.length;
    if (productLength < minProduct || productLength > maxProduct) continue;
    const f = pickCandidate(region, targetTm, minLen, maxLen, { minTm, maxTm, ...salt });
    const r = pickCandidate(revcomp(region), targetTm, minLen, maxLen, { minTm, maxTm, ...salt });
    if (!f || !r) { checks.push(`no primer met the Tm window for product ${productLength} bp`); continue; }
    const fSeq = (forwardTail + f.seq).toUpperCase();
    const rSeq = (reverseTail + r.seq).toUpperCase();
    const dimerDg = pairDimer(fSeq, rSeq).dg;
    const tmDiff = Math.abs(f.tm - r.tm);
    if (tmDiff > maxTmDiff) { checks.push(`Tm difference ${tmDiff.toFixed(1)} °C exceeds ${maxTmDiff}`); continue; }
    let extraSites = 0;
    if (specificity) {
      const fSites = findBindingSites(seq, circular, fSeq, minAnneal, {});
      const rSites = findBindingSites(seq, circular, rSeq, minAnneal, {});
      extraSites = Math.max(0, fSites.filter((s) => s.strand === 1).length - 1) + Math.max(0, rSites.filter((s) => s.strand === -1).length - 1);
    }
    const dimerPenalty = dimerDg < -5 ? (-5 - dimerDg) * 2 : 0;
    const score = Math.abs(f.tm - targetTm) + Math.abs(r.tm - targetTm) + 4 * tmDiff + f.penalty + r.penalty + 2 * extraSites + dimerPenalty;
    candidates.push({
      forward: { seq: fSeq, annealing: f.seq, tail: forwardTail, tm: f.tm, gc: gcContent(f.seq), dir: 1, properties: f.props },
      reverse: { seq: rSeq, annealing: r.seq, tail: reverseTail, tm: r.tm, gc: gcContent(r.seq), dir: -1, properties: r.props },
      productLength,
      start: mod(start0 + leftTrim, L),
      end: circular ? mod(start0 + leftTrim + region.length, L) : start0 + leftTrim + region.length,
      tmDiff,
      dimerDg,
      extraSites,
      score,
    });
  }
  candidates.sort((a, b) => a.score - b.score);
  const top = candidates.slice(0, maxCandidates);
  if (!top.length) {
    throw new Error(`No primer pair satisfied the constraints. ${checks.slice(0, 3).join('; ')}`);
  }
  const best = top[0];
  const warnings = [];
  if (best.extraSites) warnings.push(`${best.extraSites} extra binding site(s) found on the template; check specificity experimentally.`);
  if (best.dimerDg < -5) warnings.push(`Primer-dimer ΔG37 is ${best.dimerDg.toFixed(1)} kcal/mol.`);
  for (const side of ['forward', 'reverse']) {
    const hp = best[side].properties.hairpin;
    if (hp.dg < -3) warnings.push(`${side} primer may form a ${hp.stem} bp hairpin (ΔG37 ${hp.dg.toFixed(1)} kcal/mol).`);
  }
  return { ...best, candidates: top, warnings, fwdTail: forwardTail, revTail: reverseTail };
}

// ---- PCR --------------------------------------------------------------------

function pairSites(doc, fwdSeq, revSeq, opts) {
  const L = doc.seq.length;
  const fSites = findBindingSites(doc.seq, doc.circular, fwdSeq, opts.minAnneal, opts);
  const rSites = findBindingSites(doc.seq, doc.circular, revSeq, opts.minAnneal, opts);
  const products = [];
  const pairs = [[fwdSeq, revSeq, fSites, rSites], [revSeq, fwdSeq, rSites, fSites]];
  for (const [aSeq, bSeq, aS, bS] of pairs) {
    for (const a of aS.filter((s) => s.strand === 1)) {
      for (const b of bS.filter((s) => s.strand === -1)) {
        const inner = doc.circular ? mod(b.start - a.end, L) || L : b.start - a.end;
        if (inner < 0 || (!doc.circular && b.start < a.end)) continue;
        products.push({ a, b, aSeq, bSeq, inner });
      }
    }
  }
  products.sort((x, y) => x.inner - y.inner);
  return { fSites, rSites, products };
}

// Explain why PCR fails, without throwing.
export function pcrDiagnostics(doc, fwdSeq, revSeq, opts = {}) {
  const options = { minAnneal: 15, maxMismatches: 0, threePrimeExact: 1, ...opts };
  const { fSites, rSites, products } = pairSites(doc, fwdSeq, revSeq, options);
  const reasons = [];
  if (!fSites.length) reasons.push('Forward primer has no binding site (check its 3′ end and the template strand).');
  if (!rSites.length) reasons.push('Reverse primer has no binding site (check its 3′ end and the template strand).');
  if (fSites.length && rSites.length && !products.length) {
    reasons.push('Primers bind but do not face each other, or the reverse site lies upstream of the forward site.');
  }
  return { ok: products.length > 0, reasons, sites: { forward: fSites, reverse: rSites }, products };
}

// Simulate PCR on a template document with two primers (sequences 5'->3').
export function simulatePCR(doc, fwdSeq, revSeq, {
  minAnneal = 15, maxMismatches = 0, threePrimeExact = 1, maxProducts = 1, allowMismatch = true,
} = {}) {
  const opts = { minAnneal, maxMismatches: allowMismatch ? maxMismatches : 0, threePrimeExact };
  const diag = pcrDiagnostics(doc, fwdSeq, revSeq, opts);
  if (!diag.ok) {
    const e = /** @type {any} */ (new Error(`No product. ${diag.reasons.join(' ')}`));
    e.diagnostics = diag;
    throw e;
  }
  const products = diag.products.slice(0, Math.max(1, maxProducts));
  const built = products.map((p) => buildProduct(doc, p, opts));
  const warnings = [];
  if (diag.products.length > 1) warnings.push(`${diag.products.length} possible products; ${built.length === 1 ? 'the shortest was used' : `the first ${built.length} are shown`}.`);
  const mismatches = built.reduce((n, b) => n + b.mismatches, 0);
  if (mismatches) warnings.push(`The primers anneal with ${mismatches} mismatch${mismatches === 1 ? '' : 'es'}.`);
  const main = built[0].doc;
  main.history = [...(doc.history ?? [])];
  addHistory(main, 'PCR', `PCR from ${doc.name} (${built[0].start + 1}–${built[0].end}), ${main.seq.length} bp`);
  main.notes.description = `PCR product from ${doc.name}`;
  return {
    doc: main,
    products: built.map((b) => ({ doc: b.doc, start: b.start, end: b.end, length: b.doc.seq.length, mismatches: b.mismatches })),
    junctions: built.map((b) => ({ left: b.start, right: b.end })),
    warnings,
    diagnostics: diag,
  };
}

function buildProduct(doc, pair, opts) {
  const { a, b, aSeq, bSeq } = pair;
  const L = doc.seq.length;
  const inner = doc.circular ? subseq(doc.seq, a.end, b.start, true) : doc.seq.slice(a.end, b.start);
  const leftTop = a.annealed; // primer 5'->3' matches the top strand
  const rightTop = revcomp(b.annealed);
  const product = (a.tail + leftTop + inner + rightTop + revcomp(b.tail)).toUpperCase();
  const start = doc.circular ? mod(a.end - a.annealed.length - a.tail.length, L) : a.end - a.annealed.length - a.tail.length;
  const out = newDocument({ name: `${doc.name} PCR product`, seq: product, circular: false });
  const regionStart = doc.circular ? mod(a.start, L) : a.start;
  out.features = featuresInRegion(doc, regionStart, doc.circular ? mod(b.end - regionStart, L) || L : b.end - regionStart, a.tail.length);
  const fwdStart = a.tail.length;
  const revStart = product.length - b.tail.length - b.annealed.length;
  out.features.push(primerFeature('forward primer', aSeq, fwdStart, fwdStart + a.annealed.length, 1));
  out.features.push(primerFeature('reverse primer', bSeq, revStart, revStart + b.annealed.length, -1));
  out.primers = [
    { id: newId(), name: 'forward primer', seq: aSeq.toUpperCase(), description: '', extra: {} },
    { id: newId(), name: 'reverse primer', seq: bSeq.toUpperCase(), description: '', extra: {} },
  ];
  return { doc: out, start: regionStart, end: doc.circular ? mod(b.end, L) || L : b.end, mismatches: a.mismatches + b.mismatches, pair };
}

function primerFeature(name, seq, start, end, strand) {
  return {
    id: newId(), name, type: 'primer_bind', strand, color: '#5b8fd6',
    segments: [{ start, end, color: '#5b8fd6' }],
    qualifiers: [{ name: 'label', values: [{ kind: 'text', value: name }] }],
    extra: { primer: seq.toUpperCase() },
  };
}

// Add primer_bind features to a document for every primer binding site.
export function primerBindingFeatures(doc, primer, { minAnneal = 15, name = null, ...opts } = {}) {
  const seq = typeof primer === 'string' ? primer : primer.seq;
  const label = name ?? (typeof primer === 'string' ? 'primer' : primer.name);
  return findBindingSites(doc.seq, doc.circular, seq, minAnneal, opts).map((s) => primerFeature(
    label, seq, s.start, s.end, s.strand,
  ));
}

// ---- gel --------------------------------------------------------------------

// Fragments or plain sizes -> a drawable gel description.
export function gelView(fragments, { ladder = '1 kb ladder', agarose = 1, label = true } = {}) {
  const sizes = fragments.map((f) => (typeof f === 'number' ? f : f.len ?? f.length ?? 0)).filter((n) => n > 0);
  const bands = sizes
    .map((bp, index) => ({ index, bp, migration: migration(bp, agarose), label: label ? `${bp} bp` : '' }))
    .sort((a, b) => a.bp - b.bp);
  const ladderSizes = LADDERS[ladder] ?? LADDERS['1 kb ladder'];
  return {
    agarose,
    ladder,
    bands,
    ladderBands: ladderSizes.map((bp) => ({ bp, migration: migration(bp, agarose), label: `${bp}` })),
    maxBp: Math.max(0, ...sizes),
    minBp: Math.min(Infinity, ...sizes.filter((n) => n > 0)),
    list: bands.map((b) => b.bp),
  };
}

export { duplexThermo, hairpin, dimer, clampScore, LADDERS, migration };
