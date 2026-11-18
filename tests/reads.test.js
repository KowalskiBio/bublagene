import test from 'node:test';
import assert from 'node:assert/strict';
import { JSDOM } from 'jsdom';
import { placeRead } from '../js/readalign.js';
import { revcomp } from '../js/seq.js';
import { makeAB1 } from './ab1-fixture.js';

// a reproducible pseudo-random reference
function reference(n, seed = 7) {
  let s = seed;
  let out = '';
  for (let i = 0; i < n; i++) { s = (s * 1103515245 + 12345) & 0x7fffffff; out += 'ACGT'[(s >> 16) & 3]; }
  return out;
}
const REF = reference(4000);

test('a forward read is placed where it came from, with its mismatch shown', () => {
  const piece = REF.slice(1200, 1700).split('');
  piece[250] = piece[250] === 'A' ? 'C' : 'A';
  const read = placeRead(REF, { bases: piece.join(''), name: 'F' });
  assert.equal(read.reverse, false);
  assert.equal(read.start, 1200);
  assert.equal(read.chars.length, 500);
  assert.equal(read.mismatches, 1);
  assert.equal(read.map[0], 1200);
  assert.equal(read.map[499], 1699);
});

test('a reverse read is flipped onto the forward strand and its bases still map', () => {
  const piece = REF.slice(2500, 3000);
  const bases = revcomp(piece);
  const read = placeRead(REF, { bases, name: 'R' });
  assert.equal(read.reverse, true);
  assert.equal(read.start, 2500);
  assert.equal(read.chars.join(''), piece);
  assert.equal(read.map[0], 2999);        // the first called base is the last base on the forward strand
  assert.equal(read.map[499], 2500);
});

test('a deletion and an insertion in the read are carried into the row', () => {
  const piece = REF.slice(500, 1000);
  const edited = piece.slice(0, 200) + piece.slice(203, 350) + 'GG' + piece.slice(350);
  const read = placeRead(REF, { bases: edited, name: 'indel' });
  assert.equal(read.start, 500);
  assert.equal(read.chars.filter((c) => c === '-').length, 3);
  assert.equal([...read.insertions.values()].join(''), 'GG');
  assert.equal(read.map.length, edited.length);
});

test('a read from somewhere else does not match', () => {
  assert.equal(placeRead(REF, { bases: reference(500, 99), name: 'other' }), null);
});

test('opening an .ab1 with a reference open adds a read; without one it opens its own document', async () => {
  const dom = new JSDOM('<!doctype html><html><body><div id="tabs"></div><button id="openBtn"></button><button id="newBtn"></button><select id="saveAs"><option value="">Save as…</option></select><button id="undoBtn"></button><button id="redoBtn"></button><button id="themeBtn"></button><input type="file" id="fileInput" multiple><div id="viewSeg"><button data-view="map">Map</button><button data-view="seq">Sequence</button><button data-view="trace" hidden>Trace</button></div><div id="viewOpts"></div><div id="editOps"></div><div id="view"></div><div id="status"></div><div id="panelTabs"></div><div id="panel"></div><div id="toast" hidden></div></body></html>', { url: 'http://localhost/', pretendToBeVisual: true });
  const { window } = dom;
  // jsdom has no canvas; the chromatogram only needs a context that accepts drawing calls
  window.HTMLCanvasElement.prototype.getContext = () => new Proxy({}, { get: (t, k) => (k === 'measureText' ? () => ({ width: 6 }) : () => {}), set: () => true });
  Object.assign(globalThis, {
    window, document: window.document, HTMLElement: window.HTMLElement, Node: window.Node, Blob: window.Blob,
    getComputedStyle: window.getComputedStyle.bind(window),
    addEventListener: window.addEventListener.bind(window), removeEventListener: window.removeEventListener.bind(window),
    requestAnimationFrame: (cb) => setTimeout(cb, 0), cancelAnimationFrame: (id) => clearTimeout(id),
  });
  Object.defineProperty(globalThis, 'navigator', { value: window.navigator, configurable: true });
  Object.defineProperty(globalThis, 'localStorage', { value: window.localStorage, configurable: true });
  const tick = (ms = 30) => new Promise((r) => setTimeout(r, ms));
  await import(`../js/ui/app.js?reads=${Date.now()}`);
  const { openFiles } = await import('../js/ui/files.js');
  const { alignments, traces, state } = await import('../js/ui/core.js');

  const piece = REF.slice(1000, 1400);
  const ab1 = makeAB1({ bases: revcomp(piece), qualities: Array(400).fill(40), name: 'read-R' });
  // the trace comes first in the batch; the reference must still be opened before it
  const gb = `LOCUS       ref   ${REF.length} bp    DNA     linear   UNK 01-JAN-2000\nORIGIN\n${REF.match(/.{1,60}/g).map((l, i) => `${String(i * 60 + 1).padStart(9)} ${l.toLowerCase().match(/.{1,10}/g).join(' ')}`).join('\n')}\n//\n`;
  await openFiles([{ name: 'read.abi', data: ab1 }, { name: 'ref.gb', text: async () => gb }]);
  await tick(60);
  const doc = state.docs[state.active];
  assert.equal(doc.name, 'ref');
  assert.equal(state.docs.length, 1, 'the read did not become its own tab');
  const reads = alignments.get(doc.id);
  assert.equal(reads.length, 1);
  assert.equal(reads[0].reverse, true);
  assert.equal(reads[0].start, 1000);
  assert.equal(traces.get(doc.id).length, 1);
  assert.equal(traces.get(doc.id)[0].refMap[0], 1399);
  assert.ok(window.document.querySelector('#viewSeg [data-view="trace"]:not([hidden])'), 'trace view is offered');

  // a read that matches nothing opens on its own, as a chromatogram
  const stray = makeAB1({ bases: reference(300, 5), qualities: Array(300).fill(40), name: 'stray' });
  await openFiles([{ name: 'stray.ab1', data: stray }]);
  await tick(60);
  assert.equal(state.docs.length, 2);
  assert.equal(state.view, 'trace');
});

test('aligning a document that is a trace to a reference opens the Alignment view with its chromatogram', async () => {
  const { openFiles } = await import('../js/ui/files.js');
  const { alignments, traces, state } = await import('../js/ui/core.js');
  const { api } = await import('../js/ui/api/index.js');
  const { newDocument } = await import('../js/model.js');
  const { openDoc } = await import('../js/ui/documents.js');
  const tick = (ms = 30) => new Promise((r) => setTimeout(r, ms));

  // a second reference, then a read that was opened on its own (no reference was active when it came in)
  const ref = newDocument({ name: 'ref2', seq: REF, circular: false });
  const piece = REF.slice(2000, 2400);
  openDoc(newDocument({ name: 'blank', seq: '', circular: false }));
  await openFiles([{ name: 'solo.ab1', data: makeAB1({ bases: revcomp(piece), qualities: Array(400).fill(40), name: 'solo' }) }]);
  await tick();
  openDoc(ref);
  const solo = state.docs.find((d) => d.name === 'solo');
  assert.ok(solo, 'the read opened as its own document');
  state.active = state.docs.indexOf(solo);
  state.alignTarget = ref.id;
  await api.runAlignment();
  await tick();

  assert.equal(state.docs[state.active].id, ref.id, 'the reference is now the active tab');
  assert.equal(state.view, 'align');
  const rows = alignments.get(ref.id);
  assert.equal(rows.length, 1);
  assert.equal(rows[0].reverse, true);
  assert.equal(rows[0].start, 2000);
  assert.equal(traces.get(ref.id)[0].refMap[0], 2399);
});

test('a trace written the way real instruments do is read: inline small values, separate channels, unusual tag names', async () => {
  const { readAB1 } = await import('../js/formats/ab1.js');
  const { makeStandardAB1 } = await import('./ab1-fixture.js');
  const bases = 'GATTACAGATTACACCGGTT';
  const t = readAB1(makeStandardAB1({ bases, name: 'SEQ1_H10', order: 'GATC' }));
  assert.equal(t.bases, bases);
  assert.equal(t.name, 'SEQ1_H10');
  assert.equal(t.positions.length, bases.length);
  assert.equal(t.positions[0], 10 - 1 + 0, 'peak positions are 1-based in the file');
  // each called base has its own channel peak at its position: this is what the FWO_ order decides
  for (const [i, b] of [...bases].entries()) assert.equal(t.samples[b][t.positions[i] + 1] > 0 || t.samples[b][t.positions[i]] > 0, true, `peak for base ${i}`);
  assert.equal(t.tags.includes('Dye#1'), true);
});
