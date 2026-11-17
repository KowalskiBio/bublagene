// Feature clipboard, overlays (ORFs and primer sites) and typing at the caret.

import test from 'node:test';
import assert from 'node:assert/strict';
import { JSDOM } from 'jsdom';
import { newDocument, newId } from '../js/model.js';
import { copyFeatures, pasteFeatures, FEATURE_CLIP_KIND } from '../js/fragment.js';
import { orfsFor, orfFeatures, primerOverlays, isOrf, isPrimer, isOverlay } from '../js/ui/overlays.js';

const feat = (name, start, end, extra = {}) => ({
  id: newId(), name, type: 'misc_feature', strand: 1, color: '#888888',
  segments: [{ start, end }], qualifiers: [], extra: {}, ...extra,
});
const seq = (n) => 'ACGT'.repeat(Math.ceil(n / 4)).slice(0, n);

test('copied features paste at a new position, keeping their spacing', () => {
  const doc = newDocument({ seq: seq(200), circular: false });
  doc.features.push(feat('a', 10, 20), feat('b', 30, 45));
  const clip = copyFeatures(doc, doc.features.map((f) => f.id));
  assert.equal(clip.kind, FEATURE_CLIP_KIND);
  assert.equal(clip.origin, 10);
  const res = pasteFeatures(doc, 100, clip);
  assert.deepEqual(res, { added: 2, skipped: 0 });
  const [, , a2, b2] = doc.features;
  assert.deepEqual(a2.segments[0], { start: 100, end: 110 });
  assert.deepEqual(b2.segments[0], { start: 120, end: 135 });
  assert.equal(new Set(doc.features.map((f) => f.id)).size, 4, 'pasted features get new ids');
  assert.equal(doc.seq.length, 200, 'the bases are untouched');
});

test('a pasted feature that would run past the end of a linear sequence is skipped', () => {
  const doc = newDocument({ seq: seq(100), circular: false });
  doc.features.push(feat('a', 0, 10), feat('b', 50, 90));
  const clip = copyFeatures(doc, doc.features.map((f) => f.id));
  const res = pasteFeatures(doc, 40, clip);
  assert.equal(res.added, 1);
  assert.equal(res.skipped, 1);
});

test('features wrap around the origin of a circular sequence', () => {
  const doc = newDocument({ seq: seq(100), circular: true });
  doc.features.push(feat('a', 0, 20));
  const clip = copyFeatures(doc, [doc.features[0].id]);
  pasteFeatures(doc, 90, clip);
  const pasted = doc.features[1].segments[0];
  assert.equal(pasted.start, 90);
  assert.equal(pasted.end, 10, 'the end wraps to base 10');
});

test('a feature that already wraps keeps its length when copied', () => {
  const doc = newDocument({ seq: seq(100), circular: true });
  doc.features.push(feat('w', 90, 10));
  const clip = copyFeatures(doc, [doc.features[0].id]);
  assert.equal(clip.features[0].segments[0].length, 20);
  pasteFeatures(doc, 40, clip);
  assert.deepEqual(doc.features[1].segments[0], { start: 40, end: 60 });
});

test('multi-segment features keep each segment and its extras', () => {
  const doc = newDocument({ seq: seq(300), circular: false });
  doc.features.push({ ...feat('cds', 10, 20), segments: [{ start: 10, end: 20 }, { start: 50, end: 60, color: '#123456' }] });
  const clip = copyFeatures(doc, [doc.features[0].id]);
  pasteFeatures(doc, 100, clip);
  const [s1, s2] = doc.features[1].segments;
  assert.deepEqual([s1.start, s1.end, s2.start, s2.end], [100, 110, 140, 150]);
  assert.equal(s2.color, '#123456');
});

test('nothing to copy or paste is reported, not guessed', () => {
  const doc = newDocument({ seq: seq(50), circular: false });
  assert.equal(copyFeatures(doc, []), null);
  assert.throws(() => pasteFeatures(doc, 0, null), /No features/);
  assert.throws(() => pasteFeatures(newDocument({ seq: '' }), 0, { kind: FEATURE_CLIP_KIND, features: [{ segments: [] }] }), /no sequence/);
});

test('ORFs become hollow overlay bars, and the scan is cached until the sequence changes', () => {
  const orf = `ATG${'GCT'.repeat(120)}TAA`;
  const doc = newDocument({ seq: `${'C'.repeat(30)}${orf}${'C'.repeat(30)}`, circular: false });
  const orfs = orfsFor(doc, 100);
  assert.equal(orfs.length, 1);
  assert.equal(orfsFor(doc, 100), orfs, 'same answer, no rescan');
  const [bar] = orfFeatures(orfs);
  assert.ok(isOrf(bar) && isOverlay(bar) && !isPrimer(bar));
  assert.equal(bar.segments[0].start, 30);
  assert.equal(bar.name, 'ORF 121 aa');
  assert.notEqual(orfsFor({ ...doc, seq: `${doc.seq}A` }, 100), orfs, 'an edited sequence is rescanned');
});

test('every place a primer anneals becomes an overlay bar', () => {
  const site = 'GATTACAGGCTTAACCGGTA';
  const doc = newDocument({ seq: `${seq(40)}${site}${seq(40)}${site}${seq(20)}`, circular: false });
  doc.primers.push({ id: 'p1', name: 'fwd', seq: site, description: '', extra: {} });
  const bars = primerOverlays(doc);
  assert.equal(bars.length, 2);
  assert.ok(bars.every((b) => isPrimer(b) && b.primerId === 'p1' && b.type === 'primer_bind'));
  assert.notEqual(bars[0].id, bars[1].id);
});

// ---- typing in the sequence view ---------------------------------------------

test('typing, Backspace, arrows and paste edit the sequence at the caret', async () => {
  const dom = new JSDOM(`<!doctype html><html><body>
    <div id="tabs"></div><button id="openBtn"></button><button id="newBtn"></button>
    <select id="saveAs"><option value="">Save as…</option></select>
    <button id="undoBtn"></button><button id="redoBtn"></button><button id="themeBtn"></button>
    <input type="file" id="fileInput" multiple>
    <div id="viewSeg"><button data-view="map" aria-pressed="true">Map</button><button data-view="seq">Sequence</button><button data-view="trace" hidden>Trace</button></div>
    <div id="viewOpts"></div><div id="editOps"></div><div id="view"></div><div id="status"></div>
    <div id="panelTabs"></div><div id="panel"></div><div id="toast" hidden></div></body></html>`, { url: 'http://localhost/', pretendToBeVisual: true });
  const { window } = dom;
  Object.assign(globalThis, {
    window, document: window.document, HTMLElement: window.HTMLElement, Node: window.Node, Blob: window.Blob,
    addEventListener: window.addEventListener.bind(window), removeEventListener: window.removeEventListener.bind(window),
    requestAnimationFrame: (cb) => setTimeout(cb, 0), cancelAnimationFrame: (id) => clearTimeout(id),
  });
  Object.defineProperty(globalThis, 'navigator', { value: window.navigator, configurable: true });
  Object.defineProperty(globalThis, 'localStorage', { value: window.localStorage, configurable: true });
  await import(`../js/ui/app.js?editing=${Date.now()}`);
  const tick = (ms = 30) => new Promise((r) => setTimeout(r, ms));
  const $ = (s) => window.document.querySelector(s);
  const status = () => $('#status').textContent;
  const key = (k, opts = {}) => window.document.dispatchEvent(new window.KeyboardEvent('keydown', { key: k, bubbles: true, cancelable: true, ...opts }));
  await tick();

  $('#newBtn').click();
  await tick();
  $('#pasteSeq').value = '>x\nACGTACGTAC';
  [...window.document.querySelectorAll('#view button')].find((b) => b.textContent === 'Add sequence').click();
  await tick();
  assert.match(status(), /10 bp/);

  $('#viewSeg [data-view="seq"]').click();
  await tick(60);
  key('g');
  await tick();
  assert.match(status(), /10 bp/, 'with no caret, typing does nothing');

  // a click without a drag leaves the caret there (jsdom has no layout: line 0, base 3)
  const wrap = $('.seqwrap');
  const at = { bubbles: true, clientX: 62 + 8.42 * 3, clientY: 1, button: 0 };
  wrap.dispatchEvent(new window.MouseEvent('pointerdown', at));
  wrap.dispatchEvent(new window.MouseEvent('pointerup', at));
  await tick();
  assert.ok($('.caret'), 'a caret is drawn');

  key('G');
  key('G');
  await tick();
  assert.match(status(), /12 bp/);
  assert.match($('.bases').textContent, /^ACGGGTACGT$/);

  key('Backspace');
  await tick();
  assert.match(status(), /11 bp/);

  key('ArrowRight', { shiftKey: true });
  key('ArrowRight', { shiftKey: true });
  await tick();
  assert.match(status(), /Selection 5–6 · 2 bp/, 'Shift+Arrow grows a selection from the caret');

  key('T');
  await tick();
  assert.match(status(), /10 bp/, 'typing over a selection replaces it');

  const paste = new window.Event('paste', { bubbles: true, cancelable: true });
  paste.clipboardData = { getData: () => '>frag\nTT TT' };
  window.document.dispatchEvent(paste);
  await tick();
  assert.match(status(), /14 bp/, 'pasted bases land at the caret, headers and spaces dropped');

  $('#undoBtn').click();
  await tick();
  assert.match(status(), /10 bp/);

  // every panel and cloning method renders for a plain circular sequence
  const tabs = [...window.document.querySelectorAll('#panelTabs button')].map((b) => b.textContent);
  for (const name of ['Features', 'Detect', 'Mutate', 'Cloning', 'File']) assert.ok(tabs.includes(name), `${name} tab`);
  for (const name of ['Detect', 'Mutate', 'Enzymes', 'Cloning']) {
    [...window.document.querySelectorAll('#panelTabs button')].find((b) => b.textContent === name).click();
    await tick();
    assert.ok($('#panel').textContent.trim().length > 20, `${name} panel has content`);
  }
  for (const mode of ['goldengate', 'gibson', 'infusion', 'recombination', 'topo', 'restriction']) {
    const pick = window.document.getElementById('cloneMode');
    pick.value = mode;
    pick.dispatchEvent(new window.Event('change', { bubbles: true }));
    await tick();
    assert.ok($('#panel').textContent.trim().length > 40, `${mode} section has content`);
  }
  [...window.document.querySelectorAll('#panelTabs button')].find((b) => b.textContent === 'Detect').click();
  await tick();
  [...window.document.querySelectorAll('#panel button')].find((b) => b.textContent === 'Detect features').click();
  await tick();
  assert.match($('#panel').textContent, /new element|Nothing new/);
});
