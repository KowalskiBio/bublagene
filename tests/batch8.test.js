import test from 'node:test';
import assert from 'node:assert/strict';
import {
  fold, dotBracketToPairs, pairsToDotBracket, structureView, structureToText,
  structureToSvg, FOLD_MODEL, hairpinDG, duplexDG,
} from '../js/fold.js';

test('documented hairpins fold to their expected dot-bracket structures', () => {
  const a = fold('GGGGAAAACCCC');
  assert.equal(a.structure, '((((....))))');
  assert.equal(a.energy, -4.3);
  assert.equal(a.pairs.length, 4);
  const b = fold('GCGCGAAAACGCGC');
  assert.equal(b.structure, '(((((....)))))');
  const c = fold('GGGAAACCC');
  assert.equal(c.structure, '(((...)))');
  assert.equal(fold('AAAAAAAAAA').structure, '..........');
  assert.equal(fold('AAAAAAAAAA').energy, 0);
  assert.equal(a.model, FOLD_MODEL);
  assert.match(a.limits, /no pseudoknots/i);
});

test('folding is deterministic and respects the length bound', () => {
  const one = fold('GGGGAAAACCCC');
  const two = fold('GGGGAAAACCCC');
  assert.deepEqual(one, two);
  assert.throws(() => fold('A'.repeat(300), { maxLength: 250 }), /limited to 250/);
  assert.equal(fold('GGGGAAAACCCC', { maxLength: 12 }).structure, '((((....))))');
});

test('dot-bracket and pairs round-trip', () => {
  const structure = '((((....))))';
  const pairs = dotBracketToPairs(structure);
  assert.deepEqual(pairs, [[0, 11], [1, 10], [2, 9], [3, 8]]);
  assert.equal(pairsToDotBracket(12, pairs), structure);
  assert.throws(() => dotBracketToPairs('((())'), /Unbalanced/);
});

test('structure text and SVG exports describe the fold', () => {
  const seq = 'GGGGAAAACCCC';
  const { structure } = fold(seq);
  const rows = structureView(seq, structure, { width: 6 });
  assert.equal(rows.length, 2);
  assert.equal(rows[0].sequence, 'GGGGAA');
  const text = structureToText(seq, structure);
  assert.ok(text.includes('((((....))))'));
  const svg = structureToSvg(seq, structure);
  assert.ok(svg.startsWith('<svg'));
  assert.ok(svg.includes('<path'));
  assert.ok(svg.includes('>G</text>'));
});

test('duplex and hairpin ΔG are shared with the primer thermodynamics', () => {
  assert.ok(hairpinDG('GCGCGCAAAATGCGCGC').dg < 0);
  assert.ok(duplexDG('ATGCGTAC', 'GTACGCAT').dg < 0);
});
