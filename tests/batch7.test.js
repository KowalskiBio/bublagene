import test from 'node:test';
import assert from 'node:assert/strict';
import {
  findGuides, findOffTargets, guidesToFeatures, nhejEdit, hdrEdit, frameEffectReport,
  editPrimers, PAM_PRESETS,
} from '../js/crispr.js';
import { newDocument, newId } from '../js/model.js';
import { revcomp } from '../js/seq.js';

const GUIDE = 'GACCTGACCTGACCTGACCT'; // 20 nt, 60% GC

test('SpCas9 guides are found at the expected locus on the top strand', () => {
  const seq = 'TTT' + GUIDE + 'TGG' + 'AAA';
  const guides = findGuides(seq, { bothStrands: false });
  assert.equal(guides.length, 1);
  assert.equal(guides[0].sequence, GUIDE);
  assert.equal(guides[0].strand, 1);
  assert.equal(guides[0].pam, 'TGG');
  assert.equal(guides[0].start, 3);
  assert.equal(guides[0].end, 23);
  assert.equal(guides[0].cut, 20);
  assert.ok(guides[0].gc >= 20 && guides[0].gc <= 80);
});

test('both strands are searched and the guide sequence is 5′→3′', () => {
  const top = 'TTT' + GUIDE + 'TGG' + 'AAA';
  const seq = revcomp(top);
  const guides = findGuides(seq);
  const minus = guides.find((g) => g.strand === -1);
  assert.ok(minus, 'expected a guide on the minus strand');
  assert.equal(minus.sequence, GUIDE);
  assert.equal(findGuides(seq, { bothStrands: false }).length, 1);
});

test('GC bounds and configurable PAMs are respected', () => {
  const atRich = 'TTT' + 'ATATATATATATATATATAT' + 'TGG';
  assert.equal(findGuides(atRich, { minGC: 40 }).length, 0);
  const sa = 'TTT' + GUIDE + 'A' + 'AAGAGT' + 'AAA';
  const found = findGuides(sa, { enzyme: 'SaCas9' });
  assert.ok(found.some((g) => g.pam === 'AAGAGT'));
  const cas12 = 'TTTA' + GUIDE + 'AAAA';
  const five = findGuides(cas12, { enzyme: 'Cas12a' });
  assert.ok(five.some((g) => g.strand === 1 && g.sequence === GUIDE));
  assert.equal(PAM_PRESETS.SpCas9.guideLength, 20);
});

test('off-target search is explicitly heuristic and weights the seed', () => {
  const off = 'GACCTGACCTGACCTGACCA' + 'TGG'; // one mismatch at the 3′ end
  const seq = 'TTT' + GUIDE + 'TGG' + 'A'.repeat(30) + off + 'CCC';
  const hits = findOffTargets(GUIDE, seq, { maxMismatches: 2, seedLength: 8 });
  assert.ok(hits.length >= 1);
  const hit = hits.find((h) => h.sequence.startsWith(off.slice(0, 20)));
  assert.ok(hit);
  assert.equal(hit.heuristic, true);
  assert.match(hit.note, /not a cleavage prediction/i);
  assert.ok(hit.mismatches >= 1);
  assert.ok(hit.score >= hit.mismatches);
});

test('guides become crispr_site features', () => {
  const doc = newDocument({ seq: 'TTT' + GUIDE + 'TGG' });
  const guides = findGuides(doc.seq);
  doc.features = guidesToFeatures(doc, guides);
  assert.equal(doc.features.length, 1);
  assert.equal(doc.features[0].type, 'crispr_site');
  assert.equal(doc.features[0].qualifiers.find((q) => q.name === 'guide_sequence').values[0].value, GUIDE);
});

test('NHEJ edits are deterministic and report frame effects', () => {
  const seq = 'ATG' + 'GCT'.repeat(60) + 'TAA';
  const doc = newDocument({ seq });
  doc.features = [{
    id: newId(), name: 'cds', type: 'CDS', strand: 1, color: '#888',
    segments: [{ start: 0, end: seq.length }],
    qualifiers: [{ name: 'translation', values: [{ kind: 'text', value: '' }] }], extra: {},
  }];
  const guide = { start: 60, end: 80, cut: 70 };
  const a = nhejEdit(doc, guide, { seed: 7, deletion: 1 });
  const b = nhejEdit(doc, guide, { seed: 7, deletion: 1 });
  assert.equal(a.doc.seq, b.doc.seq);
  assert.equal(a.change.deleted.length, 1);
  assert.equal(a.doc.seq.length, seq.length - 1);
  assert.ok(a.frameEffects.some((f) => f.name === 'cds' && f.effect === 'frameshift'));
  const c = nhejEdit(doc, guide, { seed: 99, deletion: 1, insertion: 'AA' });
  assert.equal(c.doc.seq.length, seq.length + 1);
  assert.ok(c.doc.features.some((f) => f.name === 'NHEJ junction'));
});

test('HDR inserts a donor and edit primers flank the cut', () => {
  const seq = 'G'.repeat(50) + 'TTT' + GUIDE + 'TGG' + 'C'.repeat(50);
  const doc = newDocument({ seq });
  const guides = findGuides(doc.seq);
  const result = hdrEdit(doc, guides[0], { donor: 'ACGTACGTAC' });
  assert.ok(result.doc.seq.includes('ACGTACGTAC'));
  assert.equal(result.change.inserted, 'ACGTACGTAC');
  assert.equal(result.doc.seq.length, seq.length + 10);
  const primers = editPrimers(doc, guides[0].cut, { arm: 40 });
  assert.ok(doc.seq.includes(primers.forward.sequence));
  assert.ok(doc.seq.includes(revcomp(primers.reverse.sequence)));
});

test('frame report distinguishes in-frame and premature-stop changes', () => {
  const seq = 'ATG' + 'GCT'.repeat(30) + 'TAA';
  const before = newDocument({ seq });
  before.features = [{ id: newId(), name: 'cds', type: 'CDS', strand: 1, color: '#888', segments: [{ start: 0, end: seq.length }], qualifiers: [], extra: {} }];
  const after = newDocument({ seq: seq.slice(0, 30) + 'TAA' + seq.slice(33) });
  after.features = [{ ...before.features[0], id: newId(), segments: [{ start: 0, end: after.seq.length }] }];
  const report = frameEffectReport(before, after);
  assert.equal(report.length, 1);
  assert.ok(['frameshift', 'premature-stop', 'in-frame-change'].includes(report[0].effect));
});
