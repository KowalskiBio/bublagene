import test from 'node:test';
import assert from 'node:assert/strict';
import { primerProperties, duplexThermo, degenerateTm, expandIUPAC, hairpin, dimer, tm } from '../js/thermo.js';
import {
  findBindingSites, designPrimers, designPrimerPairs, simulatePCR, pcrDiagnostics,
  primerBindingFeatures, gelView, migration, LADDERS,
} from '../js/primers.js';
import { newDocument, featureSequence, newId } from '../js/model.js';
import { revcomp } from '../js/seq.js';

const templateDoc = (seq, circular = false) => newDocument({ seq, circular });

test('Tm fixtures agree with the SantaLucia unified model within 2 °C', () => {
  const conditions = { na: 0.05, mg: 0, oligo: 250e-9 };
  const fixtures = [
    ['GTAAAACGACGGCCAGT', 51.6],
    ['CAGGAAACAGCTATGAC', 46.0],
    ['AGCGGATAACAATTTCACACAGGA', 56.7],
    ['GATTACAGATTACACCGGTTAAGGCC', 57.4],
  ];
  for (const [seq, expected] of fixtures) {
    const value = tm(seq, conditions);
    assert.ok(Math.abs(value - expected) <= 2, `${seq}: ${value.toFixed(2)} vs ${expected}`);
  }
  assert.ok(tm('ATGCATGCATGC', { na: 0.5 }) > tm('ATGCATGCATGC', { na: 0.05 }));
  assert.ok(tm('ATGCATGCATGC', { na: 0.05, mg: 0.003 }) > tm('ATGCATGCATGC', { na: 0.05 }));
});

test('primer properties: GC, 3′ clamp, degeneracy, secondary structure', () => {
  const p = primerProperties('ATGCGCGCGCGCGCGCGCAT');
  assert.ok(p.gc > 60);
  assert.equal(p.clamp.threePrimeGc, false);
  assert.equal(p.clamp.gcInLast5, 3);
  assert.equal(primerProperties('ATGCGCGCGC').clamp.threePrimeGc, true);
  const deg = primerProperties('ATGCGNNT');
  assert.equal(deg.degeneracy, 16);
  assert.ok(deg.tmRange.max > deg.tmRange.min);
  assert.equal(expandIUPAC('ARY').length, 4);
  assert.ok(hairpin('GCGCGCAAAATGCGCGC').dg < -5);
  assert.ok(dimer('ATGCGTAC', 'GTACGCAT').dg < -5);
  assert.equal(dimer('AAAAAAAA', 'AAAAAAAA').dg, 0);
  assert.equal(duplexThermo('ATGC').length, 4);
  const d = degenerateTm('ATGCGNNT');
  assert.equal(d.count, 16);
});

test('design respects product bounds, Tm window, tails and restriction sites', () => {
  const region = 'ATGGCCTACGTTAGCTAGCATCGATCGTAGCTAGCTACGATCGTAGCATCGAT';
  const flanks = 'T'.repeat(50);
  const seq = flanks + region + 'A'.repeat(40);
  const design = designPrimerPairs(seq, false, 50, 50 + region.length, {
    targetTm: 60, minTm: 45, maxTm: 72, minProduct: 40, maxProduct: 80,
  });
  assert.equal(design.productLength, region.length);
  assert.ok(design.forward.seq.length >= 18);
  assert.ok(seq.includes(design.forward.annealing));
  assert.ok(seq.includes(revcomp(design.reverse.annealing)));
  assert.ok(Math.abs(design.forward.tm - design.reverse.tm) <= 4);

  const tailed = designPrimerPairs(seq, false, 50, 50 + region.length, {
    targetTm: 60, minTm: 45, maxTm: 72, fwdTail: 'GGGG', revTail: 'CCCC',
  });
  assert.equal(tailed.forward.seq.slice(0, 4), 'GGGG');
  assert.equal(tailed.reverse.seq.slice(0, 4), 'CCCC');
  assert.equal(tailed.productLength, region.length + 8);

  const withSite = designPrimerPairs(seq, false, 50, 50 + region.length, {
    targetTm: 60, minTm: 45, maxTm: 72, fwdRestriction: 'EcoRI',
  });
  assert.ok(withSite.forward.seq.includes('GAATTC'));
});

function randomSeq(n, seed) {
  let x = seed >>> 0;
  let out = '';
  const bases = 'ACGT';
  for (let i = 0; i < n; i++) {
    x = (1103515245 * x + 12345) >>> 0;
    out += bases[(x >>> 16) & 3];
  }
  return out;
}

test('binding sites: exact, mismatch-tolerant and strict 3′ end', () => {
  const template = randomSeq(400, 11);
  const primer = template.slice(150, 171); // 21 nt
  const exact = findBindingSites(template, false, primer, 21);
  assert.equal(exact.filter((s) => s.strand === 1).length, 1);
  assert.equal(exact.find((s) => s.strand === 1).mismatches, 0);
  assert.equal(exact.find((s) => s.strand === 1).annealed, primer);

  const mismatchPrimer = primer.slice(0, 3) + (primer[3] === 'A' ? 'C' : 'A') + primer.slice(4);
  assert.equal(findBindingSites(template, false, mismatchPrimer, 21).filter((s) => s.strand === 1).length, 0);
  const tolerated = findBindingSites(template, false, mismatchPrimer, 21, { maxMismatches: 1 }).filter((s) => s.strand === 1);
  assert.equal(tolerated.length, 1);
  assert.equal(tolerated[0].mismatches, 1);

  const threePrime = primer.slice(0, -1) + (primer.at(-1) === 'A' ? 'C' : 'A');
  assert.equal(
    findBindingSites(template, false, threePrime, 15, { maxMismatches: 3, threePrimeExact: 3 }).filter((s) => s.strand === 1).length,
    0,
  );

  const sites = primerBindingFeatures(newDocument({ seq: template }), primer);
  assert.equal(sites.length, 1);
  assert.equal(sites[0].type, 'primer_bind');
  assert.equal(sites[0].strand, 1);
});

test('PCR builds the expected product and carries primer_bind features', () => {
  const left = 'GATTACAGATTACACCGGTTAAGGCC';
  const middle = 'T'.repeat(30) + 'GAATTC' + 'A'.repeat(30);
  const right = 'CCTTGGAACCTTGGCATCGATCGTAG';
  const template = templateDoc('T'.repeat(40) + left + middle + right + 'T'.repeat(40));
  template.features = [{ id: newId(), name: 'site', type: 'misc_feature', strand: 1, color: '#888',
    segments: [{ start: 40 + left.length + 30, end: 40 + left.length + 36 }], qualifiers: [], extra: {} }];
  const fwd = 'GGATCC' + left;
  const rev = revcomp(right);
  const result = simulatePCR(template, fwd, rev);
  assert.equal(result.doc.seq, 'GGATCC' + left + middle + right);
  const binds = result.doc.features.filter((f) => f.type === 'primer_bind');
  assert.equal(binds.length, 2);
  assert.deepEqual(binds.map((f) => f.strand), [1, -1]);
  assert.equal(binds[0].name, 'forward primer');
  assert.equal(featureSequence(result.doc, result.doc.features.find((f) => f.type === 'misc_feature')), 'GAATTC');
  assert.equal(result.junctions.length, 1);
});

test('3′-mismatch primer yields no product; diagnostics explain why', () => {
  const seq = randomSeq(400, 23);
  const template = templateDoc(seq);
  const fwd = seq.slice(100, 122);
  const rev = revcomp(seq.slice(260, 282));
  const ok = simulatePCR(template, fwd, rev, { minAnneal: 18 });
  assert.ok(ok.doc.seq.length > 150);
  const three = fwd.slice(0, -1) + (fwd.at(-1) === 'A' ? 'C' : 'A');
  assert.throws(() => simulatePCR(template, three, rev, { minAnneal: 18, maxMismatches: 3, threePrimeExact: 3 }), /No product|No binding site/);
  const diag = pcrDiagnostics(template, three, rev, { minAnneal: 18, maxMismatches: 3, threePrimeExact: 3 });
  assert.equal(diag.ok, false);
  assert.ok(diag.reasons.length >= 1);
});

test('PCR across a circular origin still works', () => {
  const seq = 'ACGTTGCAAGGCTTAACCGG' + 'A'.repeat(200) + 'TTGGCCAAGGTTCCAAGGTT';
  const doc = templateDoc(seq, true);
  const fwd = seq.slice(-20, -2);
  const rev = revcomp(seq.slice(2, 20));
  const { doc: product } = simulatePCR(doc, fwd, rev);
  assert.equal(product.seq, seq.slice(-20, -2) + seq.slice(-2) + seq.slice(0, 20));
  assert.equal(product.circular, false);
});

test('gel view: sizes, ordering and ladder', () => {
  const gel = gelView([500, 1500, 100], { ladder: '1 kb ladder' });
  assert.deepEqual(gel.list, [100, 500, 1500]);
  assert.ok(gel.bands[0].migration > gel.bands[2].migration);
  assert.equal(gel.ladderBands.length, LADDERS['1 kb ladder'].length);
  assert.equal(gel.minBp, 100);
  assert.equal(gel.maxBp, 1500);
  assert.ok(migration(1000, 1) < migration(200, 1));
});
