// Nearest-neighbour thermodynamics for DNA oligonucleotides.
//
// The stack parameters are the SantaLucia (1998) unified set (ΔH kcal/mol,
// ΔS cal/mol·K). Tm uses the standard two-state expression with a monovalent
// salt correction; Mg2+ is folded into an equivalent Na+ concentration with the
// common von Ahsen approximation. Secondary-structure values are a bounded
// heuristic (contiguous stems only) intended for ranking primers, not for
// predicting folding equilibria.

import { revcomp, IUPAC } from './seq.js';

export const NN = {
  AA: [-7.9, -22.2], TT: [-7.9, -22.2],
  AT: [-7.2, -20.4], TA: [-7.2, -21.3],
  CA: [-8.5, -22.7], TG: [-8.5, -22.7],
  GT: [-8.4, -22.4], AC: [-8.4, -22.4],
  CT: [-7.8, -21.0], AG: [-7.8, -21.0],
  GA: [-8.2, -22.2], TC: [-8.2, -22.2],
  CG: [-10.6, -27.2],
  GC: [-9.8, -24.4],
  GG: [-8.0, -19.9], CC: [-8.0, -19.9],
};

const R = 1.987; // cal/mol·K
const INIT_H = 0.2;
const INIT_S = -5.7;

export function equivalentSodium({ na = 0.05, mg = 0 } = {}) {
  return Math.max(1e-4, na + 4 * Math.sqrt(Math.max(0, mg)));
}

// ΔH (kcal/mol), ΔS (cal/mol·K), ΔG37 (kcal/mol) and Tm (°C) for a duplex.
export function duplexThermo(seq, { na = 0.05, mg = 0, oligo = 250e-9, saltCorrection = true } = {}) {
  const q = String(seq).toUpperCase().replace(/[^ACGT]/g, '');
  if (q.length < 2) return { seq: q, dH: 0, dS: 0, dG37: 0, tm: 0, length: q.length };
  let dH = 0;
  let dS = 0;
  for (let i = 0; i < q.length - 1; i++) {
    const p = NN[q.slice(i, i + 2)];
    if (!p) continue;
    dH += p[0];
    dS += p[1];
  }
  dH += INIT_H;
  dS += INIT_S;
  // SantaLucia unified terminal penalties: only A/T ends are penalised.
  for (const end of [q[0], q[q.length - 1]]) {
    if (end === 'A' || end === 'T') { dH += 2.2; dS += 6.9; }
  }
  if (saltCorrection) dS += 0.368 * (q.length - 1) * Math.log(equivalentSodium({ na, mg }));
  const selfComp = q === revcomp(q);
  if (selfComp) dS += -1.4; // symmetry correction
  const ct = selfComp ? oligo : oligo / 4;
  const tm = (dH * 1000) / (dS + R * Math.log(ct)) - 273.15;
  const dG37 = dH - (310.15 * dS) / 1000;
  return { seq: q, dH, dS, dG37, tm, length: q.length, selfComplementary: selfComp };
}

export function tm(seq, options) {
  return duplexThermo(seq, options).tm;
}

// Expand IUPAC ambiguity codes. Throws above `limit` so a pathological primer
// cannot allocate without bound.
export function expandIUPAC(seq, limit = 512) {
  let out = [''];
  for (const ch of String(seq).toUpperCase()) {
    const set = IUPAC[ch];
    if (!set) throw new Error(`Not a nucleotide code: ${ch}`);
    if (set.length > 1) {
      const next = [];
      for (const prefix of out) for (const b of set) next.push(prefix + b);
      if (next.length > limit) throw new Error(`Primer degeneracy exceeds ${limit} sequences`);
      out = next;
    } else {
      out = out.map((p) => p + set);
    }
  }
  return out;
}

// Degenerate primers: Tm of every concrete sequence, reported as a range.
export function degenerateTm(seq, options = {}, limit = 512) {
  const variants = expandIUPAC(seq, limit);
  const values = variants.map((v) => tm(v, options)).sort((a, b) => a - b);
  const mean = values.reduce((n, v) => n + v, 0) / values.length;
  return { count: values.length, min: values[0], max: values[values.length - 1], mean, variants: values };
}

// ---- secondary structure (bounded heuristic) --------------------------------

function stackSum(top, bottom) {
  // top and bottom pair as top[i] vs bottom[i]; both must be canonical pairs.
  let dg = 0;
  let pairs = 0;
  for (let i = 0; i < top.length; i++) {
    const a = top[i];
    const b = bottom[i];
    if (!complementary(a, b)) return null;
    pairs++;
    if (i > 0) {
      const p = NN[top.slice(i - 1, i + 1)] ?? NN[revcomp(top.slice(i - 1, i + 1))];
      if (p) dg += p[0] - (310.15 * p[1]) / 1000;
    } else {
      dg += a === 'G' || a === 'C' ? -1.0 : -0.5;
    }
  }
  return { dg, pairs };
}

function complementary(a, b) {
  return revcomp(a) === b;
}

// Best contiguous stem-loop a single strand can form with itself. We check the
// 5' and 3' ends facing each other (the common primer hairpin) and report the
// most stable stem of 3–8 bp.
export function hairpin(seq) {
  const s = String(seq).toUpperCase().replace(/[^ACGT]/g, '');
  let best = { dg: 0, stem: 0, loop: 0, structure: '' };
  for (let loop = 3; loop <= 12; loop++) {
    for (let stem = 3; stem <= 8; stem++) {
      if (2 * stem + loop > s.length) break;
      const left = s.slice(0, stem);
      const right = s.slice(s.length - stem); // 5'->3' on the 3' arm
      const bottom = right.split('').reverse().join(''); // antiparallel partner
      if (left !== revcomp(right)) continue;
      const stack = stackSum(left, bottom);
      if (!stack) continue;
      const loopPenalty = 4.0 + 0.6 * Math.log(loop);
      const dg = stack.dg + loopPenalty;
      if (dg < best.dg) {
        best = {
          dg,
          stem,
          loop,
          structure: `${left}${'.'.repeat(loop)}${right}`,
        };
      }
    }
  }
  return best;
}

// Best contiguous duplex between two strands (a or its reverse complement vs b).
// Mismatches break a run; a run of >=4 complementary bases is reported.
export function dimer(a, b) {
  const x = String(a).toUpperCase().replace(/[^ACGT]/g, '');
  const y = String(b).toUpperCase().replace(/[^ACGT]/g, '');
  const yr = y.split('').reverse().join(''); // b aligned antiparallel to a
  let best = { dg: 0, length: 0, offset: 0, structure: '' };
  for (let offset = -(y.length - 4); offset < x.length - 3; offset++) {
    let run = '';
    let runStart = 0;
    for (let i = Math.max(0, offset); i < Math.min(x.length, offset + yr.length); i++) {
      const bBase = yr[i - offset];
      const aBase = x[i];
      if (complementary(aBase, bBase)) {
        if (!run) runStart = i;
        run += aBase;
      } else {
        if (run.length >= 4) grade(run, runStart, offset);
        run = '';
      }
    }
    if (run.length >= 4) grade(run, runStart, offset);
  }
  function grade(run, start, off) {
    // approximate: base-pair energy via stacks of the top strand
    let dg = 0;
    for (let i = 0; i < run.length; i++) {
      if (i > 0) {
        const p = NN[run.slice(i - 1, i + 1)];
        if (p) dg += p[0] - (310.15 * p[1]) / 1000;
      } else dg += -0.5;
    }
    if (dg < best.dg) {
      best = { dg, length: run.length, offset: off, structure: `${run} (${start + 1})` };
    }
  }
  return best;
}

// ---- primer properties ------------------------------------------------------

export function clampScore(seq) {
  const s = String(seq).toUpperCase().slice(-5);
  let gc = 0;
  for (const ch of s) if (ch === 'G' || ch === 'C') gc++;
  return { gcInLast5: gc, threePrimeGc: /[GC]$/.test(seq.toUpperCase()) };
}

export function primerProperties(seq, options = {}) {
  const s = String(seq).toUpperCase();
  const degenerate = /[^ACGT]/.test(s.replace(/U/g, 'T'));
  const core = s.replace(/[^ACGT]/g, '');
  const thermo = duplexThermo(core, options);
  const deg = degenerate ? degenerateTm(s, options) : null;
  const hp = hairpin(core);
  const self = dimer(core, core);
  return {
    sequence: s,
    length: s.length,
    gc: core.length ? (100 * [...core].filter((c) => c === 'G' || c === 'C').length) / core.length : 0,
    tm: thermo.tm,
    tmRange: deg ? { min: deg.min, max: deg.max, mean: deg.mean, count: deg.count } : { min: thermo.tm, max: thermo.tm, mean: thermo.tm, count: 1 },
    dH: thermo.dH,
    dS: thermo.dS,
    dG37: thermo.dG37,
    hairpin: hp,
    selfDimer: self,
    degeneracy: deg ? deg.count : 1,
    clamp: clampScore(core),
  };
}

export function pairDimer(a, b) {
  return dimer(a, b);
}
