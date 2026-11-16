import test from 'node:test';
import assert from 'node:assert/strict';
import { newDocument, newId } from '../js/model.js';
import { revcomp } from '../js/seq.js';
import { ENZYME_BY_NAME, findSites, digest } from '../js/enzymes.js';
import { readDna, writeDna } from '../js/formats/dna.js';
import {
  endsCompatible, describeEnd, pieceFromFragment, checkLigation, restrictionClone,
  assemblePieces, gibsonAssemble, cloningRecord, cloningFromJSON, cloningSidecarName, cloningSummary,
} from '../js/cloning.js';

const feat = (name, start, end, strand = 1) => ({
  id: newId(), name, type: 'misc_feature', strand, color: '#888',
  segments: [{ start, end }], qualifiers: [], extra: {},
});

const cut = (doc, names) =>
  names.flatMap((n) => findSites(doc.seq, doc.circular, ENZYME_BY_NAME[n]));

test('ends ligate when overhangs match, not when they are reverse complements', () => {
  assert.equal(endsCompatible({ type: '5', oh: 'AATT' }, { type: '5', oh: 'AATT' }), true);
  assert.equal(endsCompatible({ type: '5', oh: 'aatt' }, { type: '5', oh: 'AATT' }), true);
  assert.equal(endsCompatible({ type: '5', oh: 'AATT' }, { type: '5', oh: 'GATC' }), false);
  assert.equal(endsCompatible({ type: '5', oh: 'AAGC' }, { type: '5', oh: 'GCTT' }), false);
  assert.equal(endsCompatible({ type: '3', oh: 'TGCA' }, { type: '3', oh: 'TGCA' }), true);
  assert.equal(endsCompatible({ type: 'blunt', oh: '' }, { type: 'blunt', oh: '' }), true);
  assert.equal(endsCompatible({ type: 'blunt', oh: '' }, { type: '5', oh: 'AATT' }), false);
  assert.equal(endsCompatible({ type: 'natural', oh: '' }, { type: 'natural', oh: '' }), false);
  assert.equal(endsCompatible({ type: '5', oh: 'AATT' }, { type: '3', oh: 'AATT' }), false);
  assert.match(describeEnd({ type: '5', oh: 'AATT' }), /5′ overhang AATT/);
});

test('restriction cloning: EcoRI + BamHI, insert goes in backwards', () => {
  const vector = newDocument({
    name: 'vector.dna',
    seq: 'A'.repeat(40) + 'GAATTC' + 'C'.repeat(100) + 'GGATCC' + 'T'.repeat(30),
    circular: true,
  });
  const insert = newDocument({
    name: 'insert.dna',
    seq: 'G'.repeat(20) + 'GAATTC' + 'A'.repeat(60) + 'GGATCC' + 'G'.repeat(20),
    circular: false,
  });
  vector.features = [feat('vector core', 46, 146)];
  insert.features = [feat('insert core', 26, 86)];
  vector.features.push(feat('stuffer', 152, 182));

  const vecFrags = digest(vector.seq, vector.circular, cut(vector, ['EcoRI', 'BamHI']));
  const insFrags = digest(insert.seq, insert.circular, cut(insert, ['EcoRI', 'BamHI']));
  const backbone = [...vecFrags].sort((a, b) => b.len - a.len)[0];
  const insertPiece = insFrags[1];

  const check = checkLigation({
    backbone: { doc: vector, fragment: backbone },
    insert: { doc: insert, fragment: insertPiece },
    circular: true,
  });
  assert.deepEqual(check.orientations, ['reverse']);
  assert.equal(check.forward.compatible, false);
  assert.equal(check.forward.productSize, 172);
  assert.match(check.forward.problems.join(' '), /AATT.*GATC/);

  const { doc, orientation, junctions } = restrictionClone({
    backbone: { doc: vector, fragment: backbone },
    insert: { doc: insert, fragment: insertPiece },
    circular: true,
    name: 'construct',
  });

  const expected = 'AATTC' + 'C'.repeat(100) + 'G' + 'GATCC' + 'T'.repeat(60) + 'G';
  assert.equal(orientation, 'reverse');
  assert.equal(doc.seq, expected);
  assert.equal(doc.seq.length, 172);
  assert.equal(doc.circular, true);

  const byName = Object.fromEntries(doc.features.map((f) => [f.name, f]));
  assert.deepEqual(Object.keys(byName).sort(), ['insert core', 'vector core']);
  assert.deepEqual([byName['vector core'].segments[0].start, byName['vector core'].segments[0].end], [5, 105]);
  assert.deepEqual([byName['insert core'].segments[0].start, byName['insert core'].segments[0].end], [111, 171]);
  assert.equal(byName['insert core'].strand, -1);

  assert.equal(junctions.length, 2);
  assert.equal(junctions[0].position, 106);
  assert.ok(junctions[0].regeneratedSites.some((s) => s.name === 'BamHI' && s.pos === 106));
  assert.ok(junctions[1].regeneratedSites.some((s) => s.name === 'EcoRI' && s.pos === 0));
  assert.equal(doc.cloning.inputs[0].role, 'backbone');
  assert.equal(doc.cloning.inputs[1].orientation, 'reverse');
  assert.ok(doc.history.some((h) => h.action === 'Restriction cloning'));
  assert.ok(doc.notes.modified);
});

test('restriction cloning: a site that can ligate both ways is reported', () => {
  const vector = newDocument({
    name: 'v.dna',
    seq: 'A'.repeat(60) + 'GAATTC' + 'C'.repeat(60),
    circular: true,
  });
  const insert = newDocument({
    name: 'i.dna',
    seq: 'G'.repeat(20) + 'GAATTC' + 'T'.repeat(40) + 'GAATTC' + 'G'.repeat(20),
    circular: false,
  });
  const backbone = digest(vector.seq, true, cut(vector, ['EcoRI']))[0];
  const insertPiece = digest(insert.seq, false, cut(insert, ['EcoRI']))[1];

  const check = checkLigation({
    backbone: { doc: vector, fragment: backbone },
    insert: { doc: insert, fragment: insertPiece },
    circular: true,
  });
  assert.deepEqual(check.orientations, ['forward', 'reverse']);
  assert.equal(check.preferred, 'forward');
  assert.match(check.warnings.join(' '), /either orientation/);

  const forward = restrictionClone({
    backbone: { doc: vector, fragment: backbone },
    insert: { doc: insert, fragment: insertPiece },
    circular: true,
    orientation: 'forward',
    name: 'fwd',
  });
  const reverse = restrictionClone({
    backbone: { doc: vector, fragment: backbone },
    insert: { doc: insert, fragment: insertPiece },
    circular: true,
    orientation: 'reverse',
    name: 'rev',
  });
  assert.equal(forward.doc.seq.length, 172);
  assert.equal(reverse.doc.seq.length, 172);
  assert.notEqual(forward.doc.seq, reverse.doc.seq);
  assert.equal(
    forward.doc.seq,
    'AATTC' + 'C'.repeat(60) + 'A'.repeat(60) + 'G' + 'AATTC' + 'T'.repeat(40) + 'G',
  );
});

test('an incompatible orientation is refused with the overhangs named', () => {
  const vector = newDocument({ name: 'v', seq: 'A'.repeat(40) + 'GAATTC' + 'C'.repeat(100) + 'GGATCC' + 'T'.repeat(30), circular: true });
  const insert = newDocument({ name: 'i', seq: 'G'.repeat(20) + 'GAATTC' + 'A'.repeat(60) + 'GGATCC' + 'G'.repeat(20), circular: false });
  const backbone = [...digest(vector.seq, true, cut(vector, ['EcoRI', 'BamHI']))].sort((a, b) => b.len - a.len)[0];
  const insertPiece = digest(insert.seq, false, cut(insert, ['EcoRI', 'BamHI']))[1];
  assert.throws(
    () => restrictionClone({
      backbone: { doc: vector, fragment: backbone },
      insert: { doc: insert, fragment: insertPiece },
      circular: true,
      orientation: 'forward',
    }),
    /cannot go in forward.*AATT.*GATC/i,
  );
});

test('assemblePieces uses non-palindromic ends and rejects the wrong order', () => {
  const mk = (name, seq, leftOh, rightOh) => ({
    doc: { name, seq, circular: false, features: [], primers: [] },
    source: { name, seq, circular: false, features: [], primers: [] },
    fragment: null,
    reverse: false,
    left: { type: '5', oh: leftOh, enzyme: null },
    right: { type: '5', oh: rightOh, enzyme: null },
  });
  const a = mk('A', 'AAAA', 'AAGC', 'AATG');
  const b = mk('B', 'CCCC', 'AATG', 'AAGC');
  const { doc } = assemblePieces([a, b], { circular: true, name: 'ab' });
  assert.equal(doc.seq, 'AAAACCCC');
  assert.equal(doc.circular, true);
  assert.equal(doc.cloning.junctions.length, 2);

  const wrong = mk('C', 'GGGG', 'AAGC', 'AATG');
  assert.throws(() => assemblePieces([a, wrong], { circular: true }), /cannot join/);
});

test('Gibson assembly joins two fragments sharing terminal overlaps', () => {
  const ovA = 'ACGTACGTACGTACGT';
  const ovB = 'GGCCTTAAGGCCTTAA';
  const f1 = newDocument({ name: 'f1', seq: ovA + 'A'.repeat(30) + ovB, circular: false });
  const f2 = newDocument({ name: 'f2', seq: ovB + 'C'.repeat(30) + ovA, circular: false });
  f1.features = [feat('A block', 16, 46)];
  f2.features = [feat('C block', 16, 46)];

  const { doc, junctions, warnings } = gibsonAssemble([f1, f2], { circular: true, minOverlap: 15, name: 'gibson' });
  assert.equal(doc.seq, ovA + 'A'.repeat(30) + ovB + 'C'.repeat(30));
  assert.equal(doc.seq.length, 92);
  assert.equal(doc.circular, true);
  assert.equal(junctions.length, 2);
  assert.equal(junctions[0].position, 62);
  assert.equal(junctions[0].overhang, ovB);
  assert.equal(junctions[1].position, 92);
  assert.equal(junctions[1].overhang, ovA);
  assert.equal(warnings.length, 0);

  const byName = Object.fromEntries(doc.features.map((f) => [f.name, f]));
  assert.deepEqual([byName['A block'].segments[0].start, byName['A block'].segments[0].end], [16, 46]);
  assert.deepEqual([byName['C block'].segments[0].start, byName['C block'].segments[0].end], [62, 92]);
  assert.equal(doc.cloning.operation, 'gibson-assembly');
  assert.ok(doc.history.some((h) => h.action === 'Gibson assembly'));
});

test('Gibson assembly flips a fragment when the overlap requires it', () => {
  const ov = 'TTGACCAGTTGACCAGT';
  const f1 = newDocument({ name: 'left', seq: 'A'.repeat(20) + ov, circular: false });
  const desired = newDocument({ name: 'right', seq: ov + 'C'.repeat(20), circular: false });
  const flipped = newDocument({ name: 'right', seq: revcomp(desired.seq), circular: false });

  const { doc, inputs } = gibsonAssemble([f1, flipped], { circular: false, minOverlap: 15, name: 'linear' });
  assert.equal(doc.seq, 'A'.repeat(20) + ov + 'C'.repeat(20));
  assert.equal(doc.circular, false);
  assert.equal(inputs[1].orientation, 'reverse');
});

test('Gibson assembly reports when no overlap is found', () => {
  const f1 = newDocument({ name: 'a', seq: 'A'.repeat(40), circular: false });
  const f2 = newDocument({ name: 'b', seq: 'C'.repeat(40), circular: false });
  assert.throws(() => gibsonAssemble([f1, f2], { circular: true, minOverlap: 15 }), /No assembly found/);
});

test('cloning history is a portable JSON record', () => {
  const vector = newDocument({ name: 'pUC.dna', seq: 'A'.repeat(40) + 'GAATTC' + 'C'.repeat(100) + 'GGATCC' + 'T'.repeat(30), circular: true });
  const insert = newDocument({ name: 'gene.dna', seq: 'G'.repeat(20) + 'GAATTC' + 'A'.repeat(60) + 'GGATCC' + 'G'.repeat(20), circular: false });
  vector.features = [feat('ori', 46, 146)];
  insert.features = [feat('gene', 26, 86)];
  const backbone = [...digest(vector.seq, true, cut(vector, ['EcoRI', 'BamHI']))].sort((a, b) => b.len - a.len)[0];
  const insertPiece = digest(insert.seq, false, cut(insert, ['EcoRI', 'BamHI']))[1];
  const { doc } = restrictionClone({
    backbone: { doc: vector, fragment: backbone },
    insert: { doc: insert, fragment: insertPiece },
    circular: true,
    name: 'pUC-gene.dna',
  });

  const record = cloningRecord(doc);
  assert.equal(record.schema, 'bublagene/cloning/1');
  assert.equal(record.type, 'restriction');
  assert.equal(record.operation, 'restriction-ligation');
  assert.equal(record.product.seqLen, 172);
  assert.equal(record.inputs.length, 2);
  assert.equal(record.inputs[1].orientation, 'reverse');
  assert.equal(record.junctions.length, 2);
  assert.ok(record.junctions[0].regeneratedSites.some((s) => s.name === 'BamHI' && s.pos === 106));

  assert.deepEqual(cloningFromJSON(JSON.stringify(record)), record);
  assert.throws(() => cloningFromJSON({ schema: 'other' }), /Unsupported cloning record/);
  assert.equal(cloningSidecarName('a/construct.dna'), 'a/construct.cloning.json');
  assert.equal(cloningSidecarName('construct'), 'construct.cloning.json');

  const summary = cloningSummary(doc);
  assert.equal(summary.size, 172);
  assert.equal(summary.inputs.length, 2);
  assert.equal(summary.junctions.length, 2);
  assert.ok(summary.junctions[0].regenerated.includes('BamHI at 107'));
  assert.equal(cloningRecord(newDocument()), null);
  assert.match(doc.notes.description, /Restriction cloning: .*pUC\.dna.*gene\.dna.*reversed/);
  const round = readDna(writeDna(doc), 'round.dna');
  assert.equal(round.seq, doc.seq);
  assert.match(round.notes.description, /Restriction cloning: .*gene\.dna.*reversed/);
});

test('pieceFromFragment keeps feature coordinates and flips with the piece', () => {
  const doc = newDocument({ name: 'x', seq: 'A'.repeat(10) + 'GAATTC' + 'C'.repeat(20) + 'GGATCC' + 'T'.repeat(10), circular: false });
  doc.features = [feat('f', 20, 30)];
  const frags = digest(doc.seq, false, cut(doc, ['EcoRI', 'BamHI']));
  const middle = frags[1];
  const forward = pieceFromFragment(doc, middle);
  const reverse = pieceFromFragment(doc, middle, { reverse: true });
  assert.equal(forward.doc.seq.length, reverse.doc.seq.length);
  assert.ok(forward.doc.seq.includes('CCCCCCCCCCCCCCCCCCCCC'), 'forward piece carries the C block');
  assert.ok(revcomp(reverse.doc.seq).includes('CCCCCCCCCCCCCCCCCCCCC'), 'reverse piece carries the bottom strand');
  assert.equal(forward.left.oh, 'AATT');
  assert.equal(forward.right.oh, 'GATC');
  assert.equal(reverse.left.oh, 'GATC');
  assert.equal(reverse.right.oh, 'AATT');
  const f = forward.doc.features.find((x) => x.name === 'f');
  assert.ok(f, 'feature inside the fragment is kept');
  assert.equal(f.strand, 1);
  assert.equal(reverse.doc.features.find((x) => x.name === 'f').strand, -1);
});
