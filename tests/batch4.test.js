import test from 'node:test';
import assert from 'node:assert/strict';
import {
  homologyAssemble, goldenGateAssemble, goldenGatePieces, validateJunctions, replayAssembly,
  assemblyPlan, assemblyFromJSON, batchAssembly, designAssemblyPrimers, contentHash,
  checkReadingFrame, ASSEMBLY_SCHEMA,
} from '../js/assembly.js';
import { newDocument, newId } from '../js/model.js';
import { ENZYME_BY_NAME } from '../js/enzymes.js';

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

test('homology assembly joins three and four fragments scarlessly', () => {
  const a = randomSeq(120, 1);
  const bBody = randomSeq(100, 2);
  const cBody = randomSeq(100, 3);
  const b = a.slice(-20) + bBody;
  const c = b.slice(-20) + cBody;
  const result = homologyAssemble([a, b, c], { circular: false, minOverlap: 15, name: 'linear3' });
  assert.equal(result.doc.seq, a + bBody + cBody);
  assert.equal(result.plan.schema, ASSEMBLY_SCHEMA);
  assert.equal(result.plan.method, 'homology');
  assert.equal(result.plan.junctions.length, 2);
  assert.deepEqual(result.plan.junctions.map((j) => j.overhang.length).sort((x, y) => x - y), [20, 20]);
  const replay = replayAssembly(result.plan);
  assert.equal(replay.ok, true);
  assert.equal(replay.seq, result.doc.seq);

  const dBody = randomSeq(80, 4);
  const d = c.slice(-20) + dBody;
  const four = homologyAssemble([a, b, c, d], { circular: false });
  assert.equal(four.doc.seq, a + bBody + cBody + dBody);
  assert.equal(four.plan.fragments.length, 4);
});

test('circular homology assembly closes the circle', () => {
  const aBody = randomSeq(100, 5);
  const bBody = randomSeq(100, 6);
  const cBody = randomSeq(100, 7);
  const a = aBody + 'GGGGGGGGGGGGGGGGGGGG';
  const b = a.slice(-20) + bBody;
  const c = b.slice(-20) + cBody + a.slice(0, 20);
  const result = homologyAssemble([a, b, c], { circular: true });
  assert.equal(result.doc.circular, true);
  assert.equal(result.doc.seq.length, 320);
  assert.ok(result.doc.seq.startsWith(aBody));
  assert.equal(replayAssembly(result.plan).ok, true);
});

test('Golden Gate explicit overhangs: circular, linear and mismatch', () => {
  const f1 = { name: 'f1', seq: 'ATGAAACCC', left: 'AATG', right: 'GGTA' };
  const f2 = { name: 'f2', seq: 'TTTGGGAAA', left: 'GGTA', right: 'AATG' };
  const circular = goldenGateAssemble([f1, f2], { circular: true });
  assert.equal(circular.doc.seq, 'AATGATGAAACCC' + 'GGTATTTGGGAAA');
  assert.equal(circular.junctions.length, 2);
  assert.equal(circular.doc.assembly.method, 'golden-gate');
  assert.equal(replayAssembly(circular.plan).ok, true);

  const linear = goldenGateAssemble([f1, f2], { circular: false });
  assert.equal(linear.doc.seq, 'AATGATGAAACCC' + 'GGTATTTGGGAAA' + 'AATG');
  assert.equal(linear.junctions.length, 1);

  assert.throws(() => goldenGateAssemble([f1, { ...f2, left: 'CCCC' }], { circular: true }), /Overhangs do not match/);
});

test('Golden Gate overhangs derived from Type IIS cut sites', () => {
  const seq = 'GGTCTC' + 'A' + 'CGTT' + 'ATGCATGCAT' + 'GGGGG' + 'GAGACC' + 'TTTT';
  const piece = goldenGatePieces(seq, { enzyme: 'BsaI' });
  assert.equal(piece.seq, 'ATGCATGCAT');
  assert.equal(piece.left, 'CGTT');
  assert.equal(piece.right, 'GGGG');
  assert.equal(piece.evidence.left.type, '5');
  const f2 = { name: 'f2', seq: 'CCCCAAAA', left: 'GGGG', right: 'CGTT' };
  const result = goldenGateAssemble([piece, f2], { circular: true, enzyme: 'BsaI' });
  assert.equal(result.doc.seq, 'CGTT' + 'ATGCATGCAT' + 'GGGG' + 'CCCCAAAA');
  assert.equal(replayAssembly(result.plan).ok, true);
});

test('Golden Gate carries source features into the product', () => {
  const doc = newDocument({ name: 'src', seq: 'GGTCTC' + 'A' + 'CGTT' + 'ATGCATGCAT' + 'GGGGG' + 'GAGACC' + 'TTTT' });
  doc.features = [{ id: newId(), name: 'body bit', type: 'misc_feature', strand: 1, color: '#888',
    segments: [{ start: 15, end: 18 }], qualifiers: [], extra: {} }];
  const f2 = { name: 'f2', seq: 'CCCCAAAA', left: 'GGGG', right: 'CGTT' };
  const result = goldenGateAssemble([doc, f2], { circular: true, enzyme: 'BsaI' });
  assert.equal(result.doc.features.length, 1);
  const seq = result.doc.seq.slice(result.doc.features[0].segments[0].start, result.doc.features[0].segments[0].end);
  assert.equal(seq, 'ATG');
});

test('junction validation: forbidden sites, duplicate overhangs and frame', () => {
  const f1 = { name: 'f1', seq: 'ATGAAACCCTTT', left: 'AATG', right: 'GGTA' };
  const f2 = { name: 'f2', seq: 'GGTCTCAAAA', left: 'GGTA', right: 'AATG' };
  const result = goldenGateAssemble([f1, f2], { circular: true, enzyme: 'BsaI' });
  assert.ok(result.checks.some((c) => c.code === 'forbidden-site' && c.severity === 'error'));

  const duplicated = goldenGateAssemble([
    { name: 'a', seq: 'AAAA', left: 'GGGG', right: 'GGGG' },
    { name: 'b', seq: 'CCCC', left: 'GGGG', right: 'GGGG' },
  ], { circular: true });
  assert.ok(duplicated.checks.some((c) => c.code === 'duplicate-overhang' && c.severity === 'warning'));

  const doc = newDocument({ seq: 'ATGAAACCCGGG', circular: true });
  doc.features = [{ id: newId(), name: 'cds', type: 'CDS', strand: 1, color: '#888',
    segments: [{ start: 0, end: 11 }], qualifiers: [], extra: {} }];
  assert.equal(checkReadingFrame(doc, [{ position: 5 }]).code, 'frame-shift');
  assert.equal(checkReadingFrame(newDocument({ seq: 'ATGAAACCC', circular: true }), []), null);
});

test('plans convert to and from JSON and survive a replay', () => {
  const a = randomSeq(120, 9);
  const b = a.slice(-20) + randomSeq(90, 10);
  const result = homologyAssemble([a, b], { circular: false, name: 'plan-test' });
  const json = JSON.parse(JSON.stringify(assemblyPlan(result.doc)));
  const plan = assemblyFromJSON(json);
  assert.equal(plan.schema, ASSEMBLY_SCHEMA);
  assert.equal(plan.product.digest, contentHash(result.doc.seq));
  assert.equal(replayAssembly(plan).ok, true);
  assert.throws(() => assemblyFromJSON({ schema: 'other/1' }), /schema/);
});

test('batch assembly builds several constructs from a library', () => {
  const a = randomSeq(120, 11);
  const bBody = randomSeq(100, 12);
  const b = a.slice(-20) + bBody;
  const library = { a: { name: 'a', seq: a }, b: { name: 'b', seq: b } };
  const batch = batchAssembly(library, [
    { name: 'one', fragments: ['a', 'b'], method: 'homology', circular: false },
    { name: 'two', fragments: ['a', 'b'], method: 'homology', circular: false },
    { name: 'broken', fragments: ['a', 'missing'], method: 'homology' },
  ]);
  assert.equal(batch.ok, 2);
  assert.equal(batch.failed, 1);
  assert.ok(batch.results.find((r) => r.name === 'broken').error.includes('Unknown fragment'));
});

test('assembly primers add the neighbouring overlap tails', () => {
  const a = randomSeq(150, 13);
  const b = randomSeq(150, 14);
  const primers = designAssemblyPrimers([a, b], { overlap: 20, targetTm: 60 });
  assert.equal(primers.length, 2);
  assert.equal(primers[0].forward.tail, b.slice(-20));
  assert.equal(primers[0].reverse.tail, b.slice(0, 20));
  assert.ok(primers[0].forward.sequence.startsWith(primers[0].forward.tail));
  assert.ok(primers[0].reverse.sequence.startsWith(primers[0].reverse.tail));
});
