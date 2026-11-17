// Property-based tests (fast-check) for parsers, locations, digests, PCR and
// assembly. Every generator is bounded, so runs stay fast and deterministic.

import test from 'node:test';
import assert from 'node:assert/strict';
import fc from 'fast-check';
import { parseLocation, formatLocation, readGenBank, writeGenBank } from '../js/formats/genbank.js';
import { findSites, digest, ENZYME_BY_NAME } from '../js/enzymes.js';
import { designPrimerPairs, simulatePCR } from '../js/primers.js';
import { homologyAssemble } from '../js/assembly.js';
import { newDocument, newId } from '../js/model.js';
import { revcomp, rangeLength } from '../js/seq.js';

test('location parse/format round-trips random segment sets', () => {
  fc.assert(fc.property(
    fc.integer({ min: 20, max: 200 }),
    fc.array(fc.tuple(fc.integer({ min: 0, max: 190 }), fc.integer({ min: 1, max: 8 })), { minLength: 1, maxLength: 3 }),
    fc.boolean(),
    (L, raw, negative) => {
      const segs = [];
      for (const [start, len] of raw) {
        const s = Math.min(start, L - 1);
        const e = Math.min(L, s + len);
        if (e > s && !segs.some((x) => s < x.end && e > x.start)) segs.push({ start: s, end: e });
      }
      if (!segs.length) return true;
      segs.sort((a, b) => a.start - b.start);
      const f = { strand: negative ? -1 : 1, segments: segs };
      const parsed = parseLocation(formatLocation(f, L), L);
      if (!parsed) return false;
      return parsed.strand === f.strand && JSON.stringify(parsed.segments.map((x) => [x.start, x.end])) === JSON.stringify(segs.map((x) => [x.start, x.end]));
    },
  ), { numRuns: 80 });
});

test('digest fragments always add up to the sequence length', () => {
  const seqArb = fc.array(fc.constantFrom('A', 'C', 'G', 'T'), { minLength: 20, maxLength: 400 }).map((a) => a.join(''));
  fc.assert(fc.property(seqArb, fc.array(fc.constantFrom('EcoRI', 'BamHI', 'HindIII', 'SmaI', 'BsaI'), { maxLength: 5 }), (seq, names) => {
    const cuts = names.flatMap((n) => findSites(seq, false, ENZYME_BY_NAME[n]));
    const linear = digest(seq, false, cuts);
    const circular = digest(seq, true, cuts.map((c) => ({ ...c })));
    return linear.reduce((n, f) => n + f.len, 0) === seq.length && circular.reduce((n, f) => n + f.len, 0) === seq.length;
  }), { numRuns: 60 });
});

test('designed PCR products carry both primers and the target', () => {
  fc.assert(fc.property(
    fc.array(fc.constantFrom('A', 'C', 'G', 'T'), { minLength: 180, maxLength: 260 }).map((a) => a.join('')),
    (seq) => {
      const region = { start: 40, end: seq.length - 40 };
      const design = designPrimerPairs(seq, false, region.start, region.end, { targetTm: 55, minTm: 40, maxTm: 75, maxTmDiff: 20 });
      const { doc } = simulatePCR(newDocument({ seq }), design.forward.seq, design.reverse.seq, { minAnneal: 15 });
      return doc.seq.startsWith(design.forward.seq) && doc.seq.endsWith(revcomp(design.reverse.seq));
    },
  ), { numRuns: 20 });
});

test('homology assembly reconstructs random overlap chains', () => {
  fc.assert(fc.property(
    fc.array(fc.constantFrom('A', 'C', 'G', 'T'), { minLength: 40, maxLength: 80 }).map((a) => a.join('')),
    fc.array(fc.constantFrom('A', 'C', 'G', 'T'), { minLength: 40, maxLength: 80 }).map((a) => a.join('')),
    fc.array(fc.constantFrom('A', 'C', 'G', 'T'), { minLength: 40, maxLength: 80 }).map((a) => a.join('')),
    fc.integer({ min: 16, max: 30 }),
    (aBody, bBody, cBody, ov) => {
      const a = aBody;
      const b = a.slice(-Math.min(ov, a.length)) + bBody;
      const c = b.slice(-Math.min(ov, b.length)) + cBody;
      const result = homologyAssemble([a, b, c], { circular: false, minOverlap: 15 });
      const expected = a + b.slice(Math.min(ov, a.length)) + c.slice(Math.min(ov, b.length));
      return result.doc.seq === expected;
    },
  ), { numRuns: 25 });
});

test('GenBank write/read preserves random feature geometry and qualifiers', () => {
  fc.assert(fc.property(
    fc.integer({ min: 30, max: 300 }),
    fc.array(fc.tuple(fc.integer({ min: 0, max: 290 }), fc.integer({ min: 1, max: 20 })), { minLength: 1, maxLength: 6 }),
    fc.boolean(),
    (L, raw, negative) => {
      const seq = 'ATGC'.repeat(Math.ceil(L / 4)).slice(0, L);
      const doc = newDocument({ name: 'prop', seq, circular: true });
      const segs = [];
      for (const [s, len] of raw) {
        const start = Math.min(s, L - 1);
        const end = Math.min(L, start + len);
        if (end > start && !segs.some((x) => start < x.end && end > x.start)) segs.push({ start, end });
      }
      if (!segs.length) return true;
      segs.sort((a, b) => a.start - b.start);
      doc.features = [{
        id: newId(), name: 'misc_feature', type: 'misc_feature', strand: negative ? -1 : 1, color: '#888',
        segments: segs, qualifiers: [{ name: 'note', values: [{ kind: 'text', value: 'n1' }, { kind: 'text', value: 'n2' }] }], extra: {},
      }];
      const [back] = readGenBank(writeGenBank(doc));
      const f = back.features[0];
      const geometry = JSON.stringify(f.segments.map((x) => [x.start, x.end]));
      const wanted = JSON.stringify(segs.map((x) => [x.start, x.end]));
      return f.strand === (negative ? -1 : 1) && geometry === wanted && f.qualifiers.find((q) => q.name === 'note').values.length === 2;
    },
  ), { numRuns: 40 });
});
