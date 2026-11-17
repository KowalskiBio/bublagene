// Feature detection, recombination, TOPO, In-Fusion, mutagenesis and the history tree.

import test from 'node:test';
import assert from 'node:assert/strict';
import { newDocument, newId, codingSequence } from '../js/model.js';
import { revcomp, translate } from '../js/seq.js';
import { detectFeatures, hitToFeature, libraryFromDocs, normalizeLibrary, BUILTIN_LIBRARY } from '../js/detect.js';
import { recombine, topoClone, findRecombinationSites, siteKey } from '../js/recombination.js';
import { infusionAssemble, cloningSummary } from '../js/cloning.js';
import { aminoAcidChange, designMutagenesis } from '../js/mutagenesis.js';
import { buildTree } from '../js/ui/historytree.js';

const feat = (name, start, end, strand = 1, type = 'misc_feature') => ({
  id: newId(), name, type, strand, color: '#888888', segments: [{ start, end }], qualifiers: [], extra: {},
});
const pad = 'GCGC'.repeat(20);

// ---- detection ---------------------------------------------------------------

test('detection finds DNA elements on both strands and protein tags in any frame', () => {
  const flag = 'ATGGATTACAAGGATGACGACGATAAGTAA';
  const doc = newDocument({ seq: `${pad}TAATACGACTCACTATAGG${pad}${revcomp(flag)}${pad}`, circular: false });
  const hits = detectFeatures(doc);
  const by = Object.fromEntries(hits.map((h) => [h.name, h]));
  assert.equal(by['T7 promoter'].start, 80);
  assert.equal(by['T7 promoter'].strand, 1);
  assert.equal(by['FLAG tag'].strand, -1, 'the tag sits on the reverse strand');
  assert.ok(!('Enterokinase site' in by), 'a short site inside a longer tag is part of the tag');
});

test('a mismatch budget finds near matches, and none finds only exact ones', () => {
  const near = 'AATTGTGAGCGGATAACAATA'; // lac operator, last base changed
  const doc = newDocument({ seq: `${pad}${near}${pad}`, circular: false });
  assert.equal(detectFeatures(doc, BUILTIN_LIBRARY, { maxMismatchRate: 0 }).some((h) => h.name === 'lac operator'), false);
  const hit = detectFeatures(doc, BUILTIN_LIBRARY, { maxMismatchRate: 0.05 }).find((h) => h.name === 'lac operator');
  assert.equal(hit.mismatches, 1);
  assert.ok(hit.identity < 100);
});

test('hits wrap around the origin of a circular sequence', () => {
  const doc = newDocument({ seq: `ACTCACTATAGG${pad}TAATACG`, circular: true });
  const hit = detectFeatures(doc).find((h) => h.name === 'T7 promoter');
  assert.ok(hit && hit.end < hit.start, 'the end wraps');
});

test('elements already annotated are not suggested again, and accepted hits become features', () => {
  const doc = newDocument({ seq: `${pad}TAATACGACTCACTATAGG${pad}`, circular: false });
  const [hit] = detectFeatures(doc);
  const f = hitToFeature(hit);
  assert.equal(f.name, 'T7 promoter');
  assert.equal(f.type, 'promoter');
  doc.features.push(f);
  assert.equal(detectFeatures(doc).length, 0);
});

test('features of other documents teach the scanner long elements', () => {
  const element = 'ATGCGTACGTTAGCCGATCGATCGGCTAAGCTTAGGCTA';
  const source = newDocument({ name: 'known', seq: `${pad}${element}${pad}`, circular: false });
  source.features.push(feat('myGene', 80, 80 + element.length, 1, 'CDS'));
  const lib = libraryFromDocs([source]);
  assert.equal(lib.length, 1);
  const target = newDocument({ seq: `${pad.slice(4)}${revcomp(element)}${pad}`, circular: false });
  const hit = detectFeatures(target, lib)[0];
  assert.equal(hit.name, 'myGene');
  assert.equal(hit.strand, -1);
  assert.equal(hit.source, 'document');
});

test('library input is validated', () => {
  const lib = normalizeLibrary([{ name: 'ok', seq: 'acgtacgtacgtacgt' }, { name: 'short', seq: 'ACGT' }, { name: 'bad', seq: 'ACGTXXXXXXXXXXXXXX' }, { seq: 'ACGTACGTACGTACGT' }, { name: 'tag', kind: 'protein', seq: 'DYKDDDDK' }]);
  assert.deepEqual(lib.map((e) => e.name), ['ok', 'tag']);
});

// ---- recombination -----------------------------------------------------------

const core1 = 'TTGTACAAAAAAGCAGG';
const core2 = 'CTTTCTTGTACAAAG';
const entry = () => {
  const d = newDocument({ name: 'entry', circular: true, seq: `AAAAAAAAAAGGAAA${core1}CCTTTATGGGGCCCTAAGGTTT${core2}CCAAATTTTTTTTTT` });
  d.features.push(feat('attL1', 10, 37), feat('attL2', 49, 74), feat('gene', 27, 49, 1, 'CDS'));
  return d;
};
const dest = () => {
  const d = newDocument({ name: 'dest', circular: true, seq: `CCCCCCCCCCTTTTT${core1}GGAAACCAAACCAAACCAAGGG${core2}TTCCCGGGGGGGGGG` });
  d.features.push(feat('attR1', 10, 37), feat('attR2', 59, 84), feat('AmpR', 0, 8));
  return d;
};

test('site keys pair att1 with att1 whatever the arms, and lox by full name', () => {
  assert.equal(siteKey('attL1'), siteKey('attR1'));
  assert.notEqual(siteKey('attL1'), siteKey('attL2'));
  assert.equal(siteKey('loxP'), siteKey('LOXP'));
});

test('recombination swaps the stretch between two matching sites, through their shared core', () => {
  const r = recombine(entry(), dest());
  const expected = `TTCCCGGGGGGGGGGCCCCCCCCCCTTTTT${core1}CCTTTATGGGGCCCTAAGGTTT${core2}`;
  assert.equal(r.doc.seq, expected);
  assert.ok(r.doc.circular);
  const names = r.doc.features.map((f) => f.name);
  assert.ok(names.includes('AmpR'), 'the acceptor backbone keeps its features');
  assert.ok(names.includes('attB1') && names.includes('attB2'), 'an LR reaction makes attB sites');
  assert.ok(!names.includes('attR1'), 'the old sites are consumed');
  assert.equal(r.byproduct.seq.length, entry().seq.length + dest().seq.length - r.doc.seq.length);
  const ampR = r.doc.features.find((f) => f.name === 'AmpR').segments[0];
  assert.equal(r.doc.seq.slice(ampR.start, ampR.end), 'CCCCCCCC');
  assert.equal(r.doc.cloning.operation, 'site-specific-recombination');
  assert.match(cloningSummary(r.doc).junctions[0].kind, /bp recombined core/);
});

test('recombination reads attB sites by sequence, and says what is wrong when it cannot', () => {
  const fragment = newDocument({ seq: `GG${'ACAAGTTTGTACAAAAAAGCAGGCT'}ATGAAATAA${'ACCACTTTGTACAAGAAAGCTGGGT'}GG`, circular: false });
  assert.deepEqual(findRecombinationSites(fragment).map((s) => s.name), ['attB1', 'attB2']);
  assert.throws(() => recombine(entry(), newDocument({ seq: 'ACGT'.repeat(30), circular: true })), /No matching/);
  const linearDest = dest();
  linearDest.circular = false;
  assert.throws(() => recombine(entry(), linearDest), /circular/);
});

// ---- TOPO --------------------------------------------------------------------

test('TOPO-style cloning inserts at the chosen site, in either orientation', () => {
  const vector = newDocument({ name: 'vec', seq: 'ACGTACGTACGTACGTACGT', circular: true });
  vector.features.push(feat('lacZ', 0, 5), feat('AmpR', 15, 20));
  const insert = newDocument({ name: 'ins', seq: 'GATTACAGATTACA', circular: false });
  insert.features.push(feat('tag', 0, 7));
  const fwd = topoClone(vector, insert, { kind: 'TA', range: { start: 10, end: 10 } });
  assert.equal(fwd.doc.seq, `ACGTACGTAC${insert.seq}GTACGTACGT`);
  assert.equal(fwd.doc.features.find((f) => f.name === 'tag').segments[0].start, 10);
  assert.equal(fwd.doc.features.find((f) => f.name === 'AmpR').segments[0].start, 15 + insert.seq.length);
  const rev = topoClone(vector, insert, { kind: 'blunt', range: { start: 10, end: 10 }, orientation: 'reverse' });
  assert.equal(rev.doc.seq.slice(10, 24), revcomp(insert.seq));
  assert.equal(rev.doc.features.find((f) => f.name === 'tag').strand, -1);
  assert.ok(fwd.warnings.some((w) => /not directional/.test(w)));
});

test('directional TOPO needs the CACC end, and a selection is replaced', () => {
  const vector = newDocument({ name: 'vec', seq: 'ACGTACGTACGTACGTACGT', circular: true });
  const bad = newDocument({ name: 'bad', seq: 'GGGGATGAAA', circular: false });
  assert.throws(() => topoClone(vector, bad, { kind: 'directional', range: { start: 4, end: 8 } }), /CACC/);
  assert.throws(() => topoClone(vector, newDocument({ name: 'r', seq: revcomp('CACCATGAAA') }), { kind: 'directional', range: { start: 4, end: 8 } }), /other strand/);
  const good = newDocument({ name: 'good', seq: 'CACCATGAAA', circular: false });
  const r = topoClone(vector, good, { kind: 'directional', range: { start: 4, end: 8 } });
  assert.equal(r.doc.seq, `ACGT${good.seq}ACGTACGTACGT`);
  assert.throws(() => topoClone(vector, good, { kind: 'TA' }), /Select where/);
});

// ---- In-Fusion ---------------------------------------------------------------

test('In-Fusion assembles homologous ends and labels the product as such', () => {
  const ov = 'GATTACAGGCTTAACCGGTA';
  const a = newDocument({ name: 'F1', seq: `TTTTCCCCAAAAGGGGTTTTACGTTGCA${ov}`, circular: false });
  const b = newDocument({ name: 'F2', seq: `${ov}CCGGAATTCCGGTTAAGGCCTTAA`, circular: false });
  const r = infusionAssemble([{ doc: a }, { doc: b }], { circular: false });
  assert.equal(r.doc.seq, `TTTTCCCCAAAAGGGGTTTTACGTTGCA${ov}CCGGAATTCCGGTTAAGGCCTTAA`);
  assert.equal(r.doc.cloning.operation, 'in-fusion-assembly');
  assert.match(r.doc.name, /In-Fusion/);
  assert.ok(r.doc.history.some((h) => h.action === 'In-Fusion assembly'));
});

// ---- mutagenesis -------------------------------------------------------------

const orfDoc = () => {
  const doc = newDocument({ name: 'gene', seq: `${'CC'.repeat(10)}ATGGCTAGCAAAGGAGAAGAATTATAA${'GG'.repeat(10)}`, circular: true });
  doc.features.push(feat('orf', 20, 47, 1, 'CDS'));
  return doc;
};

test('an amino-acid change picks the codon that needs the fewest base changes', () => {
  const doc = orfDoc();
  const ch = aminoAcidChange(doc, doc.features[0], 2, 'E');
  assert.equal(ch.label, 'A2E');
  assert.equal(ch.oldCodon, 'GCT');
  assert.equal(ch.start, 23);
  assert.match(ch.newCodon, /^GA[AG]$/);
  assert.throws(() => aminoAcidChange(doc, doc.features[0], 2, 'A'), /already/);
  assert.throws(() => aminoAcidChange(doc, doc.features[0], 99, 'A'), /no residue/);
});

test('a reverse-strand CDS is edited on the top strand', () => {
  const orf = `ATGGCTAGCAAAGGATAA`;
  const doc = newDocument({ name: 'rev', seq: `${'CC'.repeat(10)}${revcomp(orf)}${'GG'.repeat(10)}`, circular: true });
  doc.features.push(feat('orf', 20, 38, -1, 'CDS'));
  const ch = aminoAcidChange(doc, doc.features[0], 2, 'E');
  assert.equal(ch.oldCodon, 'GCT');
  const r = designMutagenesis(doc, ch, { method: 'overlap' });
  const cds = r.mutant.features.find((f) => f.name === 'orf');
  assert.equal(translate(codingSequence(r.mutant, cds)).slice(0, 3), 'MES');
});

test('both primer designs build the same mutant, with primers that carry the change', () => {
  const doc = orfDoc();
  const ch = aminoAcidChange(doc, doc.features[0], 2, 'E');
  for (const method of ['overlap', 'backtoback']) {
    const r = designMutagenesis(doc, ch, { method });
    const cds = r.mutant.features.find((f) => f.name === 'orf');
    assert.equal(translate(codingSequence(r.mutant, cds)).slice(0, 3), 'MES');
    assert.equal(r.mutant.seq.length, doc.seq.length);
    assert.ok(r.forward.seq.includes(ch.newBases) || r.reverse.seq.includes(revcomp(ch.newBases)));
    assert.ok(r.mutant.primers.length >= 2);
  }
  const overlap = designMutagenesis(doc, ch, { method: 'overlap' });
  assert.equal(overlap.reverse.seq, revcomp(overlap.forward.seq), 'an overlapping pair is fully complementary');
  const b2b = designMutagenesis(doc, ch, { method: 'backtoback' });
  assert.ok(b2b.forward.seq.startsWith(ch.newBases), 'the forward primer starts at the change');
});

test('insertions and deletions change the length, and bad requests are refused', () => {
  const doc = orfDoc();
  const ins = designMutagenesis(doc, { start: 26, end: 26, newBases: 'GGGAAA' }, { method: 'backtoback' });
  assert.equal(ins.kind, 'insertion');
  assert.equal(ins.mutant.seq.length, doc.seq.length + 6);
  const del = designMutagenesis(doc, { start: 23, end: 26, newBases: '' }, { method: 'overlap' });
  assert.equal(del.kind, 'deletion');
  assert.equal(del.mutant.seq.length, doc.seq.length - 3);
  assert.throws(() => designMutagenesis({ ...doc, circular: false }, { start: 1, end: 2, newBases: 'A' }), /circular/);
  assert.throws(() => designMutagenesis(doc, { start: 20, end: 21, newBases: 'A' }), /same/);
  assert.throws(() => designMutagenesis(doc, { start: 5, end: 5, newBases: '' }), /Nothing/);
});

// ---- history tree ------------------------------------------------------------

test('the history tree follows open source files back through their own history', () => {
  const ov = 'GATTACAGGCTTAACCGGTA';
  const a = newDocument({ name: 'A', seq: `TTTTCCCCAAAAGGGGTTTTACGT${ov}`, circular: false });
  const b = newDocument({ name: 'B', seq: `${ov}CCGGAATTCCGGTTAAGGCC`, circular: false });
  const ab = infusionAssemble([{ doc: a }, { doc: b }], { circular: false }).doc;
  ab.name = 'AB';
  const c = newDocument({ name: 'C', seq: `${ab.seq.slice(-20)}ACGTACGTACGTACGTACGT`, circular: false });
  const abc = infusionAssemble([{ doc: ab }, { doc: c }], { circular: false }).doc;
  abc.name = 'ABC';
  const tree = buildTree(abc, [a, b, ab, c, abc]);
  assert.equal(tree.operation, 'In-Fusion assembly');
  const names = tree.children.map((n) => n.name);
  assert.deepEqual(names, ['AB', 'C']);
  assert.deepEqual(tree.children[0].children.map((n) => n.name), ['A', 'B'], 'AB is expanded into what it was made from');
  assert.equal(tree.children[0].openId, ab.id);
});
