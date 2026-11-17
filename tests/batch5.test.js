import test from 'node:test';
import assert from 'node:assert/strict';
import { alignGlobal, alignLocal, alignBanded, alignPair, dotPlot, formatAlignment } from '../js/align.js';
import { assembleContigs, overlapScore, consensusAt } from '../js/contig.js';
import { callVariants, variantsToFeatures, featuresToVariants, normalizeVariant, applyVariants, formatVariant } from '../js/variants.js';
import { newDocument } from '../js/model.js';
import { revcomp } from '../js/seq.js';

test('global alignment matches a hand-checked Needleman-Wunsch result', () => {
  const same = alignGlobal('GATTACA', 'GATTACA');
  assert.equal(same.score, 14);
  assert.equal(same.identity, 100);
  assert.equal(same.matches, 7);
  // the classic worked example: match 1, mismatch -1, linear gap -1 => score 0
  const wiki = alignGlobal('GATTACA', 'GCATGCU', { match: 1, mismatch: -1, gapOpen: 0, gapExtend: 1 });
  assert.equal(wiki.score, 0);
  assert.equal(wiki.alignedA.length, wiki.alignedB.length);
  assert.ok(formatAlignment(same).includes('100.0%'));
});

test('affine gaps, local alignment and banded alignment', () => {
  const ins = alignGlobal('ACGT', 'ACGGT');
  assert.equal(ins.matches, 4);
  assert.equal(ins.gaps, 1);
  assert.equal(ins.identity, 80);
  const del = alignGlobal('ACGGT', 'ACGT');
  assert.equal(del.gaps, 1);
  const local = alignLocal('TTGATTACAGG', 'GATTACA');
  assert.equal(local.alignedA, 'GATTACA');
  assert.equal(local.alignedB, 'GATTACA');
  assert.equal(local.score, 14);
  assert.equal(local.startA, 2);
  const banded = alignBanded('TTGATTACAGG', 'TTGATTACAGG', { band: 2 });
  assert.equal(banded.identity, 100);
  assert.equal(banded.score, 22);
});

test('dot plot has a diagonal for identical sequences', () => {
  const seq = 'ATGCGTACGTTAGCTAGCTA';
  const plot = dotPlot(seq, seq, { window: 6, step: 1, threshold: 0.9 });
  assert.ok(plot.points.length >= seq.length - 5);
  for (let i = 0; i + 6 <= seq.length; i++) assert.ok(plot.points.some(([x, y]) => x === i && y === i), `diagonal ${i}`);
  assert.equal(plot.best, 1);
});

test('contig assembly merges overlapping reads and calls consensus', () => {
  const a = 'ATGCGTACGTTAGCTAGCTACCGATCGATCGA';
  const b = a.slice(-20) + 'TTGGAACCTTGGCCAA';
  const c = b.slice(-18) + 'GGCCAATTAA';
  const contigs = assembleContigs([a, b, c], { minOverlap: 15, maxMismatchRate: 0.05 });
  assert.equal(contigs.length, 1);
  assert.equal(contigs[0].seq, a + b.slice(20) + c.slice(18));
  assert.equal(contigs[0].reads.length, 3);
  assert.ok(contigs[0].depth > 1);
  assert.ok(overlapScore(a, b, { minOverlap: 15 }).length >= 20);
});

test('quality-weighted consensus resolves a base difference', () => {
  const ref = 'ATGCGTACGT';
  const readA = { name: 'a', seq: 'ATGCGTACGT', qualities: Array(10).fill(40) };
  const readB = { name: 'b', seq: 'ATGCGTTCGT', qualities: [10, 10, 10, 10, 10, 10, 5, 10, 10, 10] };
  assert.equal(consensusAt(ref, [readA, readB], 0), 'ATGCGTACGT');
});

test('variant calling: SNPs, insertions, deletions, MNVs and normalisation', () => {
  const snp = callVariants('ATGCGTACGT', 'ATGCGTTCGT');
  assert.equal(snp.variants.length, 1);
  assert.deepEqual([snp.variants[0].type, snp.variants[0].position, snp.variants[0].ref, snp.variants[0].alt], ['SNP', 6, 'A', 'T']);
  assert.equal(formatVariant(snp.variants[0]), '7 A>T');

  const ins = callVariants('ATGCGT', 'ATGACGT');
  assert.equal(ins.variants.length, 1);
  assert.equal(ins.variants[0].type, 'INS');
  assert.equal(ins.variants[0].position, 3);
  assert.equal(ins.variants[0].alt, 'A');

  const del = callVariants('ATGACGT', 'ATGCGT');
  assert.equal(del.variants.length, 1);
  assert.equal(del.variants[0].type, 'DEL');
  assert.equal(del.variants[0].position, 3);

  const mnv = callVariants('ATGCGT', 'ATAAGT');
  assert.equal(mnv.variants.length, 1);
  assert.equal(mnv.variants[0].type, 'MNV');
  assert.equal(mnv.variants[0].ref, 'GC');
  assert.equal(mnv.variants[0].alt, 'AA');

  const normal = normalizeVariant({ type: 'DEL', position: 4, ref: 'A', alt: '' }, 'GAAAAT');
  assert.equal(normal.position, 1);
  assert.equal(applyVariants('GAAAAT', [normal]), 'GAAAT');
  assert.equal(applyVariants('ATGCGTACGT', snp.variants), 'ATGCGTTCGT');
});

test('variants round-trip through features', () => {
  const doc = newDocument({ seq: 'ATGCGTACGT' });
  const calls = callVariants(doc.seq, 'ATGCGTTCGT').variants;
  doc.features = variantsToFeatures(doc, calls);
  assert.equal(doc.features.length, 1);
  assert.equal(doc.features[0].type, 'variation');
  const back = featuresToVariants(doc);
  assert.deepEqual(back.map((v) => [v.type, v.position, v.ref, v.alt]), calls.map((v) => [v.type, v.position, v.ref, v.alt]));
});

test('circular variant calling finds a variant at the origin', () => {
  const ref = 'AAACCCGGGTTT';
  const query = 'T' + ref.slice(1);
  const result = callVariants(ref, query, { circular: true });
  assert.equal(result.circular, true);
  assert.ok(result.variants.some((v) => v.position === 0 && v.type === 'SNP' && v.ref === 'A' && v.alt === 'T'));
});

test('reverse-complemented reads assemble in the right orientation', () => {
  const a = 'ATGCGTACGTTAGCTAGCTACCGATCGATCGA';
  const b = a.slice(-20) + 'TTGGAACCTTGGCCAA';
  const contigs = assembleContigs([a, revcomp(b)], { minOverlap: 15, maxMismatchRate: 0.05 });
  assert.equal(contigs.length, 1);
  assert.equal(contigs[0].seq, a + b.slice(20));
});
