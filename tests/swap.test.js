import test from 'node:test';
import assert from 'node:assert/strict';
import { featureSequence, newDocument, newId, swapFeature } from '../js/model.js';
import { readPastedSequence, revcomp } from '../js/seq.js';

const feat = (name, start, end, strand = 1, type = 'misc_feature') => ({
  id: newId(), name, type, strand, color: '#888', segments: [{ start, end }],
  qualifiers: [{ name: 'note', values: [{ kind: 'text', value: `${name} note` }] }], extra: {},
});
const doc40 = () => {
  const doc = newDocument({ name: 'd', seq: 'A'.repeat(10) + 'C'.repeat(10) + 'G'.repeat(10) + 'T'.repeat(10), circular: true });
  const target = feat('target', 10, 20);
  doc.features.push(feat('before', 0, 5), target, feat('after', 30, 40));
  return { doc, target };
};

test('a swap puts the new bases in place and the feature grows to cover them', () => {
  const { doc, target } = doc40();
  const r = swapFeature(doc, target.id, 'ATATATATATATATATATAT');
  assert.equal(doc.seq.length, 50);
  assert.equal(doc.seq.slice(10, 30), 'ATATATATATATATATATAT');
  assert.deepEqual([r.start, r.end, r.delta, r.removed], [10, 30, 10, 0]);
  const t = doc.features.find((f) => f.name === 'target');
  assert.deepEqual(t.segments[0], { start: 10, end: 30 });
  assert.equal(t.qualifiers[0].values[0].value, 'target note', 'its annotation is kept');
});

test('features either side keep pointing at the same bases, and the order of the list is kept', () => {
  const { doc, target } = doc40();
  swapFeature(doc, target.id, 'AT');
  assert.deepEqual(doc.features.map((f) => f.name), ['before', 'target', 'after']);
  const after = doc.features.find((f) => f.name === 'after');
  assert.deepEqual(after.segments[0], { start: 22, end: 32 });
  assert.equal(featureSequence(doc, after), 'T'.repeat(10));
  assert.equal(doc.features.find((f) => f.name === 'before').segments[0].start, 0);
});

test('a feature on the reverse strand is written in as the reverse complement, so it still reads as pasted', () => {
  const doc = newDocument({ name: 'r', seq: 'AAAAAAAAAACCCCCCCCCCGGGGGGGGGG', circular: false });
  const f = feat('rev', 10, 20, -1);
  doc.features.push(f);
  swapFeature(doc, f.id, 'GATTACAGAT');
  assert.equal(doc.seq.slice(10, 20), revcomp('GATTACAGAT'));
  assert.equal(featureSequence(doc, doc.features.find((x) => x.name === 'rev')), 'GATTACAGAT');
});

test('features wholly inside are removed and counted; ones that overlap an end are trimmed', () => {
  const { doc, target } = doc40();
  doc.features.push(feat('inside', 12, 16), feat('overlaps', 5, 12));
  const r = swapFeature(doc, target.id, 'GG');
  assert.equal(r.removed, 1);
  assert.equal(doc.features.some((f) => f.name === 'inside'), false);
  const ov = doc.features.find((f) => f.name === 'overlaps');
  assert.equal(ov.segments[0].start, 5);
  assert.ok(ov.segments[0].end <= 10, 'trimmed back to the edit point');
});

test('a swap that cannot be done says why and leaves the sequence alone', () => {
  const { doc, target } = doc40();
  const before = doc.seq;
  assert.throws(() => swapFeature(doc, target.id, ''), /Nothing to swap/);
  assert.throws(() => swapFeature(doc, 'nope', 'A'), /no longer exists/);
  const joined = { ...feat('joined', 0, 3), segments: [{ start: 0, end: 3 }, { start: 6, end: 9 }] };
  doc.features.push(joined);
  assert.throws(() => swapFeature(doc, joined.id, 'A'), /one stretch/);
  const wrapping = { ...feat('wrap', 35, 5) };
  doc.features.push(wrapping);
  assert.throws(() => swapFeature(doc, wrapping.id, 'A'), /across the origin/);
  assert.equal(doc.seq, before);
});

test('a pasted sequence is read from FASTA, GenBank-style numbering and RNA', () => {
  assert.deepEqual(readPastedSequence('>some seq\natgc gtac\n   61 ttgca\n'), { seq: 'ATGCGTACTTGCA', ignored: 0 });
  assert.deepEqual(readPastedSequence('AUGGCC'), { seq: 'ATGGCC', ignored: 0 });
  assert.deepEqual(readPastedSequence('ATG-CC*'), { seq: 'ATGCC', ignored: 0 });
  assert.deepEqual(readPastedSequence('ATGEL'), { seq: 'ATG', ignored: 2 });
  assert.deepEqual(readPastedSequence('ATGRYN'), { seq: 'ATGRYN', ignored: 0 });
});

test('the sequence dialog shows the current bases, enables Apply once something is pasted, and hands back clean bases', async () => {
  const { JSDOM } = await import('jsdom');
  const dom = new JSDOM('<!doctype html><html><body></body></html>', { url: 'http://localhost/' });
  const { window } = dom;
  Object.assign(globalThis, {
    window, document: window.document, HTMLElement: window.HTMLElement, Node: window.Node,
    addEventListener: window.addEventListener.bind(window), removeEventListener: window.removeEventListener.bind(window),
  });
  const { showSequenceDialog } = await import('../js/ui/seqdialog.js');
  let got = null;
  showSequenceDialog({
    title: 'Swap “X”', current: { label: 'Current', seq: 'ACGT'.repeat(40) }, confirmLabel: 'Swap',
    notes: ({ seq }) => (seq.length % 3 ? ['frame shifts'] : []),
    onConfirm: (s) => { got = s; },
  });
  const q = (sel) => window.document.querySelector(sel);
  assert.ok(q('.seqmodal'));
  assert.match(q('.seqpreview').textContent, /1ACGTACGTAC GTACGTACGT/, 'the current sequence is shown in numbered blocks');
  assert.equal(q('.btn.primary').disabled, true, 'nothing to apply yet');
  const ta = q('.seqtext');
  ta.value = '>x\nauggc c\n12 gE';
  ta.dispatchEvent(new window.Event('input', { bubbles: true }));
  assert.equal(q('.btn.primary').disabled, false);
  assert.match(q('.seqinfo').textContent, /1 letter that are not bases|1 letters? that are not bases/);
  assert.match(q('.seqinfo').textContent, /frame shifts/);
  assert.match(q('.seqpanehead .sub').textContent, /bp/);
  q('.btn.primary').click();
  assert.equal(got, 'ATGGCCG');
  assert.equal(q('.seqmodal'), null, 'the dialog closes');
});
