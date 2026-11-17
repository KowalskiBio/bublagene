// CRISPR guide finding, heuristic off-target search and edit simulation.
//
// Off-target scoring is an explicitly labelled heuristic: PAM-proximal (seed)
// mismatches are weighted more than distal ones, and no claim is made about
// cleavage probability. Edit simulation is deterministic for a given seed.

import { revcomp, iupacRegexSource, gcContent, mod, subseq } from './seq.js';
import { newDocument, cloneDoc, newId, addHistory, replaceRange, codingSequence, translateFeature, featureSpan } from './model.js';
import { mulberry32, randomInt, randomBases } from './rng.js';
import { designPrimerPairs } from './primers.js';

export const PAM_PRESETS = {
  SpCas9: { pam: 'NGG', side: '3', guideLength: 20 },
  SaCas9: { pam: 'NNGRRT', side: '3', guideLength: 21 },
  Cas12a: { pam: 'TTTV', side: '5', guideLength: 20 },
  custom: { pam: 'NGG', side: '3', guideLength: 20 },
};

const GENE_KEYS = new Set(['SpCas9', 'SaCas9', 'Cas12a']);

function pamRegex(pam) {
  return new RegExp(`^${iupacRegexSource(pam)}$`);
}

// Find guides on both strands. The returned coordinates are 0-based,
// end-exclusive guide spans; `cut` is 3 nt from the PAM for SpCas9/SaCas9.
export function findGuides(seq, {
  enzyme = 'SpCas9', pam = null, side = null, guideLength = null,
  minGC = 20, maxGC = 80, bothStrands = true, circular = false, maxGuides = 5000,
  maxOffTargets = null, offTargetMismatches = 4,
} = {}) {
  const preset = PAM_PRESETS[enzyme] ?? PAM_PRESETS.custom;
  const pamPattern = (pam ?? preset.pam).toUpperCase();
  const pamSide = side ?? preset.side;
  const glen = guideLength ?? preset.guideLength;
  const plen = pamPattern.length;
  if (glen < 10) throw new Error('Guide length must be at least 10 nt');
  const re = pamRegex(pamPattern);
  const L = seq.length;
  const out = [];

  const scan = (hay, strand) => {
    for (let i = 0; i + glen + plen <= hay.length; i++) {
      const guide = hay.slice(i, i + glen);
      const pamSeq = pamSide === '3' ? hay.slice(i + glen, i + glen + plen) : hay.slice(i - plen, i);
      if (pamSide === '5' && i < plen) continue;
      if (!re.test(pamSeq)) continue;
      const gc = gcContent(guide);
      if (gc < minGC || gc > maxGC) continue;
      const guideStart = pamSide === '3' ? i : i;
      if (!circular && (guideStart < 0 || guideStart + glen > L)) continue;
      let start;
      let end;
      if (strand === 1) {
        start = mod(guideStart, L);
        end = mod(guideStart + glen, L) || L;
      } else {
        start = mod(L - (guideStart + glen), L);
        end = mod(L - guideStart, L) || L;
      }
      const pamStart = pamSide === '3'
        ? (strand === 1 ? mod(guideStart + glen, L) : mod(L - (guideStart + glen + plen), L))
        : (strand === 1 ? mod(guideStart - plen, L) : mod(L - guideStart, L));
      const cut = mod(strand === 1 ? guideStart + glen - 3 : L - (guideStart + glen - 3), L);
      out.push({
        sequence: guide,
        strand, start, end, pam: pamSeq, pamStart,
        cut,
        gc: Math.round(gc * 10) / 10,
        enzyme,
      });
      if (out.length >= maxGuides) return;
    }
  };

  const hay = circular ? seq + seq.slice(0, glen + plen) : seq;
  scan(hay, 1);
  if (bothStrands) scan(revcomp(circular ? seq + seq.slice(0, glen + plen) : seq), -1);
  // de-duplicate the wrap-around copy
  const seen = new Set();
  const guides = out.filter((g) => {
    const k = `${g.strand}:${g.start}:${g.end}`;
    if (seen.has(k)) return false;
    seen.add(k);
    return true;
  });
  if (maxOffTargets != null) {
    for (const g of guides) g.offTargets = findOffTargets(g.sequence, seq, { maxMismatches: offTargetMismatches, circular, pam: pamPattern, side: pamSide, guideLength: glen, limit: maxOffTargets }).filter((o) => !(o.start === g.start && o.strand === g.strand && o.mismatches === 0));
  }
  return guides.sort((a, b) => a.start - b.start || a.strand - b.strand);
}

// Heuristic off-target search. Returns PAM-adjacent near-matches with mismatch
// counts; `heuristic: true` is part of every record.
export function findOffTargets(guide, seq, {
  maxMismatches = 4, seedLength = 8, circular = false, pam = 'NGG', side = '3', guideLength = null, limit = 50,
} = {}) {
  const g = guide.toUpperCase();
  const glen = guideLength ?? g.length;
  const plen = pam.length;
  const re = pamRegex(pam);
  const L = seq.length;
  const out = [];
  const scan = (hay, strand) => {
    for (let i = 0; i + glen + plen <= hay.length; i++) {
      const pamSeq = side === '3' ? hay.slice(i + glen, i + glen + plen) : hay.slice(i - plen, i);
      if (side === '5' && i < plen) continue;
      if (!re.test(pamSeq)) continue;
      const candidate = hay.slice(i, i + glen);
      let mismatches = 0;
      let seedMismatches = 0;
      for (let k = 0; k < glen; k++) {
        if (candidate[k] !== g[k]) {
          mismatches++;
          const fromPam = side === '3' ? glen - 1 - k : k;
          if (fromPam < seedLength) seedMismatches++;
        }
      }
      if (mismatches > maxMismatches) continue;
      const start = strand === 1 ? mod(i, L) : mod(L - (i + glen), L);
      out.push({
        sequence: candidate,
        strand, start, end: mod(start + glen, L) || L,
        pam: pamSeq,
        mismatches,
        seedMismatches,
        score: mismatches + seedMismatches * 2,
        heuristic: true,
        note: 'Heuristic seed/mismatch score; not a cleavage prediction.',
      });
      if (out.length >= limit) return;
    }
  };
  const hay = circular ? seq + seq.slice(0, glen + plen) : seq;
  scan(hay, 1);
  scan(revcomp(hay), -1);
  return out.sort((a, b) => a.score - b.score || a.mismatches - b.mismatches);
}

export function guidesToFeatures(doc, guides) {
  return guides.map((g) => ({
    id: newId(), name: `guide ${g.sequence.slice(0, 8)}`, type: 'crispr_site', strand: g.strand, color: '#7a5fb0',
    segments: [{ start: g.start, end: g.end, color: '#7a5fb0' }],
    qualifiers: [
      { name: 'guide_sequence', values: [{ kind: 'text', value: g.sequence }] },
      { name: 'pam', values: [{ kind: 'text', value: g.pam }] },
      { name: 'note', values: [{ kind: 'text', value: `${g.enzyme ?? 'CRISPR'} guide; off-target scores are heuristic` }] },
    ],
    extra: { cut: g.cut },
  }));
}

// ---- edit simulation --------------------------------------------------------

function cutBounds(guide, L) {
  if (Number.isFinite(guide.cut)) return { start: guide.cut, end: guide.cut };
  return { start: guide.start, end: guide.end };
}

function addJunctionFeature(doc, at, label, note) {
  doc.features.push({
    id: newId(), name: label, type: 'misc_feature', strand: 1, color: '#c06a6a',
    segments: [{ start: Math.max(0, at - 1), end: Math.min(doc.seq.length, at + 1), color: '#c06a6a' }],
    qualifiers: [{ name: 'note', values: [{ kind: 'text', value: note }] }], extra: {},
  });
}

// Non-homologous end joining: a deterministic deletion (and optional insertion)
// at the cut. The same seed always produces the same edit.
export function nhejEdit(doc, guide, { seed = 1, deletion = null, insertion = null } = {}) {
  const rng = mulberry32(seed);
  const { start } = cutBounds(guide, doc.seq.length);
  const del = deletion ?? randomInt(rng, 1, 8);
  const ins = insertion ?? (rng() < 0.5 ? '' : randomBases(rng, randomInt(rng, 1, 3)));
  const end = Math.min(doc.seq.length, start + del);
  const product = cloneDoc(doc);
  const before = doc.seq.length;
  replaceRange(product, start, end, ins);
  addJunctionFeature(product, start, 'NHEJ junction', `NHEJ at cut ${start + 1}: deleted ${end - start} bp, inserted ${ins.length} bp (seed ${seed})`);
  addHistory(product, 'CRISPR NHEJ', `cut at ${start + 1}, deletion ${end - start} bp, insertion ${ins.length} bp`);
  product.notes.description = `NHEJ edit of ${doc.name}`;
  const change = { start, end, deleted: doc.seq.slice(start, end), inserted: ins, delta: product.seq.length - before };
  return {
    doc: product,
    change,
    frameEffects: frameEffectReport(doc, product),
    deterministic: { seed, deletion: del, insertion: ins },
    method: 'nhej',
  };
}

// Homology-directed repair: insert (or replace with) a donor sequence.
export function hdrEdit(doc, guide, { donor = '', replaceLength = 0, seed = 1 } = {}) {
  const { start } = cutBounds(guide, doc.seq.length);
  const end = Math.min(doc.seq.length, start + Math.max(0, replaceLength));
  const product = cloneDoc(doc);
  replaceRange(product, start, end, String(donor).toUpperCase());
  addJunctionFeature(product, start, 'HDR edit', `HDR at cut ${start + 1}: replaced ${end - start} bp with ${donor.length} bp donor`);
  addHistory(product, 'CRISPR HDR', `cut at ${start + 1}, donor ${donor.length} bp`);
  product.notes.description = `HDR edit of ${doc.name}`;
  return {
    doc: product,
    change: { start, end, deleted: doc.seq.slice(start, end), inserted: String(donor).toUpperCase(), delta: product.seq.length - doc.seq.length },
    frameEffects: frameEffectReport(doc, product),
    deterministic: { seed, donorLength: donor.length, replaceLength },
    method: 'hdr',
  };
}

// Compare CDS translations before and after an edit.
export function frameEffectReport(beforeDoc, afterDoc) {
  const before = new Map((beforeDoc.features ?? []).filter((f) => f.type === 'CDS').map((f) => [f.name, f]));
  const report = [];
  for (const f of afterDoc.features ?? []) {
    if (f.type !== 'CDS') continue;
    const old = before.get(f.name);
    const afterProtein = translateFeature(afterDoc, f);
    const beforeProtein = old ? translateFeature(beforeDoc, old) : null;
    const lengthChange = old ? f.segments.reduce((n, s) => n + Math.max(0, s.end - s.start), 0) - old.segments.reduce((n, s) => n + Math.max(0, s.end - s.start), 0) : null;
    let effect = 'unchanged';
    if (!old) effect = 'new';
    else if (beforeProtein !== afterProtein) {
      if (lengthChange % 3 !== 0) effect = 'frameshift';
      else if (afterProtein.slice(0, -1).includes('*') && !beforeProtein.slice(0, -1).includes('*')) effect = 'premature-stop';
      else effect = 'in-frame-change';
    }
    report.push({
      name: f.name,
      effect,
      lengthChange,
      beforeLength: beforeProtein?.length ?? null,
      afterLength: afterProtein.length,
      beforeProtein: beforeProtein ?? null,
      afterProtein,
    });
  }
  return report;
}

// PCR primers flanking the cut, for screening or overlap-extension HDR donors.
export function editPrimers(doc, position, { arm = 40, targetTm = 60, minLen = 18, maxLen = 30 } = {}) {
  const L = doc.seq.length;
  const start = Math.max(0, position - arm);
  const end = Math.min(L, position + arm);
  if (end - start < 2 * minLen) throw new Error('The region around the cut is too short for two primers');
  const design = designPrimerPairs(doc.seq, false, start, end, { targetTm, minTm: targetTm - 15, maxTm: targetTm + 30, minLen, maxLen, maxTmDiff: 20 });
  return {
    position,
    forward: { sequence: design.forward.seq, tm: design.forward.tm },
    reverse: { sequence: design.reverse.seq, tm: design.reverse.tm },
    productLength: design.productLength,
  };
}
