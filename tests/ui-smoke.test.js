// Headless smoke test for the browser wiring: loads the real app module in
// jsdom, opens documents, and renders every new panel at least once.

import test from 'node:test';
import assert from 'node:assert/strict';
import { JSDOM } from 'jsdom';

const tick = (ms = 20) => new Promise((r) => setTimeout(r, ms));

function markup() {
  return `<!doctype html><html><body>
    <div id="tabs"></div>
    <button id="openBtn"></button><button id="newBtn"></button>
    <select id="saveAs"><option value="">Save as…</option><option value="gb">gb</option></select>
    <button id="undoBtn"></button><button id="redoBtn"></button><button id="themeBtn"></button>
    <input type="file" id="fileInput" multiple>
    <div id="viewSeg"><button data-view="map" aria-pressed="true">Map</button><button data-view="seq">Sequence</button><button data-view="trace" hidden>Trace</button></div>
    <div id="viewOpts"></div><div id="editOps"></div><div id="view"></div><div id="status"></div>
    <div id="panelTabs"></div><div id="panel"></div>
    <div id="toast" hidden></div>
  </body></html>`;
}

test('the app renders maps, panels, alignment, CRISPR and structure in jsdom', async () => {
  const dom = new JSDOM(markup(), { url: 'http://localhost/', pretendToBeVisual: true });
  const { window } = dom;
  window.HTMLCanvasElement.prototype.getContext = () => null; // jsdom draws nothing; the views skip drawing without a context
  globalThis.window = window;
  globalThis.document = window.document;
  Object.defineProperty(globalThis, 'navigator', { value: window.navigator, configurable: true });
  Object.defineProperty(globalThis, 'localStorage', { value: window.localStorage, configurable: true });
  globalThis.HTMLElement = window.HTMLElement;
  globalThis.Node = window.Node;
  globalThis.Blob = window.Blob;
  globalThis.addEventListener = window.addEventListener.bind(window);
  globalThis.removeEventListener = window.removeEventListener.bind(window);
  globalThis.requestAnimationFrame = (cb) => setTimeout(cb, 0);
  globalThis.cancelAnimationFrame = (id) => clearTimeout(id);

  await import(`../js/ui/app.js?smoke=${Date.now()}`);
  await tick();
  const $ = (sel) => window.document.querySelector(sel);
  const byText = (tag, text) => [...window.document.querySelectorAll(tag)].find((e) => e.textContent.trim() === text);

  // welcome -> example document
  const openExample = byText('button', 'Open the example fragment');
  assert.ok(openExample, 'welcome button');
  openExample.click();
  await tick(40);
  assert.ok($('#view svg'), 'map rendered');
  assert.ok($('#panel').textContent.includes('Features'), 'feature panel rendered');

  // second document so alignment has a target
  ($('#view .welcome') ? byText('button', 'Open the example fragment') : null)?.click();
  if (window.document.querySelectorAll('#tabs .tab').length < 2) {
    // open the example through the same path the welcome card uses
    byText('button', 'Open the example fragment')?.click();
  }
  await tick(40);
  assert.ok(window.document.querySelectorAll('#tabs .tab').length >= 1, 'tabs');

  const panel = async (name) => { byText('#panelTabs button', name)?.click(); await tick(20); };

  // add a second (empty) document so the alignment panel has a target
  window.document.getElementById('newBtn').click();
  await tick(30);
  assert.ok(window.document.querySelectorAll('#tabs .tab').length >= 2, 'two tabs');

  // Alignment panel
  await panel('Alignment');
  const alignBtn = byText('button', 'Align');
  if (alignBtn && !alignBtn.disabled) {
    alignBtn.click();
    await tick(60);
    assert.ok($('#panel').textContent.includes('No alignment yet.') || $('#panel').textContent.includes('identity') || $('#panel').textContent.includes('second document'), 'alignment panel');
  }

  // CRISPR
  await panel('CRISPR');
  byText('button', 'Find guides')?.click();
  await tick(30);
  assert.ok($('#panel').textContent.match(/guide\(s\)|No guide/), 'crispr panel');

  // Structure (on the example document)
  window.document.querySelectorAll('#tabs .tab')[0].click();
  await tick(30);
  await panel('Structure');
  byText('button', 'Fold sequence')?.click();
  await tick(120);
  assert.ok($('#panel').textContent.includes('ΔG') || $('#panel').textContent.includes('No fold') || $('#panel').textContent.includes('limited') || $('#panel').textContent.includes('Nothing to fold'), 'structure panel');

  // Features + multi-select
  await panel('Features');
  const checks = [...window.document.querySelectorAll('#panel input[type=checkbox]')];
  assert.ok(checks.length >= 1, 'feature checkboxes');
  checks[0].click();
  await tick();
  assert.ok($('#panel').textContent.includes('Selected ('), 'multi-select section');

  // File panel shows project actions
  await panel('File');
  assert.ok($('#panel').textContent.includes('Save project'), 'project buttons');

  // Fragment copy/paste through the context menu
  window.document.querySelectorAll('#tabs .tab')[0].click();
  await tick(20);
  await panel('Features');
  const featureRow = window.document.querySelector('#panel .row');
  featureRow?.click();
  await tick(10);
  const bpBefore = Number(($('#status b')?.textContent ?? '0').replace(/,/g, ''));
  const openMenu = () => $('#view').dispatchEvent(new window.MouseEvent('contextmenu', { bubbles: true, cancelable: true, clientX: 10, clientY: 10 }));
  openMenu();
  await tick(10);
  [...window.document.querySelectorAll('.ctxmenu .item')].find((b) => b.textContent.includes('Copy fragment'))?.click();
  openMenu();
  await tick(10);
  [...window.document.querySelectorAll('.ctxmenu .item')].find((b) => b.textContent.includes('Paste fragment'))?.click();
  await tick(30);
  const bpAfter = Number(($('#status b')?.textContent ?? '0').replace(/,/g, ''));
  assert.ok(bpAfter > bpBefore, `paste increased length (${bpBefore} -> ${bpAfter})`);

  // Sequence view with the six-frame track
  const six = [...window.document.querySelectorAll('#viewOpts label')].find((l) => l.textContent.includes('Six frames'));
  six?.querySelector('input')?.click();
  await tick();
  window.document.querySelector('#viewSeg [data-view="seq"]').click();
  await tick(30);
  assert.ok($('#view .seqwrap'), 'sequence view');
  assert.ok($('#view .seqline .aa'), 'six-frame rows');

  dom.window.close();
});
