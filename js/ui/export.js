// Getting a map or plot out of the page: a standalone SVG, a PNG, a one-page
// PDF, or the sheet that is printed. Exports always use the light palette on
// white, whatever theme the page is in, so they read on paper.

import { h } from './dom.js';
import { imagePdf } from '../pdf.js';
import { circularMap, linearMap } from './map.js';
import { featureSpan } from '../model.js';

const LIGHT = {
  '--paper': '#ffffff', '--surface': '#ffffff', '--sunk': '#eceeec', '--ink': '#16211f', '--muted': '#5e6d6a',
  '--line': '#c9d2ce', '--accent': '#0d6e66', '--accent-soft': '#d7eae7', '--select': '#f2d492', '--select-line': '#b8860b',
  '--danger': '#9b2c2c', '--mismatch': '#b45309',
};

// A copy of an SVG that stands on its own: CSS variables replaced by colours,
// a white background, fonts named, and a size.
export function exportableSvg(svg) {
  const clone = /** @type {SVGSVGElement} */ (svg.cloneNode(true));
  const resolve = (v) => v.replace(/var\((--[\w-]+)\)/g, (_, name) => LIGHT[name] ?? '#000000');
  for (const node of [clone, ...clone.querySelectorAll('*')]) {
    for (const a of [...node.attributes]) if (a.value.includes('var(')) node.setAttribute(a.name, resolve(a.value));
  }
  const vb = (clone.getAttribute('viewBox') ?? '0 0 600 400').split(/\s+/).map(Number);
  const [, , w, hgt] = vb;
  clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
  clone.setAttribute('width', String(w));
  clone.setAttribute('height', String(hgt));
  clone.setAttribute('font-family', 'Helvetica, Arial, sans-serif');
  clone.removeAttribute('style');
  const bg = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
  bg.setAttribute('width', '100%');
  bg.setAttribute('height', '100%');
  bg.setAttribute('fill', '#ffffff');
  clone.insertBefore(bg, clone.firstChild);
  return { text: `<?xml version="1.0" encoding="UTF-8"?>\n${new XMLSerializer().serializeToString(clone)}`, width: w, height: hgt };
}

function drawToCanvas({ text, width, height }, scale = 2) {
  return new Promise((resolve, reject) => {
    const img = new Image();
    img.onload = () => {
      const canvas = document.createElement('canvas');
      canvas.width = Math.round(width * scale);
      canvas.height = Math.round(height * scale);
      const ctx = canvas.getContext('2d');
      ctx.fillStyle = '#ffffff';
      ctx.fillRect(0, 0, canvas.width, canvas.height);
      ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
      resolve(canvas);
    };
    img.onerror = () => reject(new Error('The browser could not draw this image'));
    img.src = `data:image/svg+xml;charset=utf-8,${encodeURIComponent(text)}`;
  });
}

const toBlob = (canvas, type, quality) => new Promise((resolve, reject) => canvas.toBlob((b) => (b ? resolve(b) : reject(new Error('Could not encode the image'))), type, quality));

export const svgBlob = (exp) => new Blob([exp.text], { type: 'image/svg+xml' });
export const pngBlob = async (exp) => toBlob(await drawToCanvas(exp, 2), 'image/png');
export async function pdfBlob(exp, { title = '', subtitle = '' } = {}) {
  const canvas = await drawToCanvas(exp, 2);
  const jpeg = new Uint8Array(await (await toBlob(canvas, 'image/jpeg', 0.92)).arrayBuffer());
  return new Blob([imagePdf(jpeg, { width: canvas.width, height: canvas.height, title, subtitle })], { type: 'application/pdf' });
}

export function download(blob, name) {
  const url = URL.createObjectURL(blob);
  const a = h('a', { href: url, download: name });
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 2000);
}

// ---- the printed sheet --------------------------------------------------------

const numbered = (seq, width = 60) => {
  const lines = [];
  for (let i = 0; i < seq.length; i += width) {
    const chunk = seq.slice(i, i + width).match(/.{1,10}/g).join(' ');
    lines.push(`${String(i + 1).padStart(7)}  ${chunk}`);
  }
  return lines.join('\n');
};

/**
 * Everything worth having on paper: the map, the features, the primers and the
 * whole sequence in numbered blocks. Built on demand, so it is complete however
 * much of the sequence is scrolled into view.
 */
export function printSheet(doc, { sites = [], style = 'auto', overlays = [] } = {}) {
  const L = doc.seq.length;
  const circular = style === 'auto' ? doc.circular : style === 'circular';
  const svg = (circular ? circularMap : linearMap)(doc, { sites, overlays });
  const exp = exportableSvg(svg);
  const sheet = h('div', { class: 'printsheet' });
  sheet.append(h('h1', {}, doc.name), h('p', {}, `${L.toLocaleString('en-US')} bp · ${doc.circular ? 'circular' : 'linear'}${doc.notes?.organism ? ` · ${doc.notes.organism}` : ''}`));
  const mapBox = h('div', { class: 'printmap' });
  mapBox.innerHTML = exp.text.replace(/^<\?xml[^>]*>\s*/, '').replace(/width="[\d.]+" height="[\d.]+"/, 'width="100%"');
  sheet.append(mapBox);
  if (doc.features.length) {
    const t = h('table', {}, h('tr', {}, ['Feature', 'Type', 'Position', 'Strand'].map((x) => h('th', {}, x))));
    for (const f of [...doc.features].sort((a, b) => featureSpan(a).start - featureSpan(b).start)) {
      const span = featureSpan(f);
      t.append(h('tr', {}, h('td', {}, f.name), h('td', {}, f.type), h('td', {}, `${span.start + 1}–${span.end || L}`), h('td', {}, f.strand === -1 ? '−' : f.strand === 1 ? '+' : '')));
    }
    sheet.append(h('h2', {}, 'Features'), t);
  }
  if (doc.primers?.length) {
    const t = h('table', {}, h('tr', {}, ['Primer', 'Sequence (5′→3′)', 'Length'].map((x) => h('th', {}, x))));
    for (const p of doc.primers) t.append(h('tr', {}, h('td', {}, p.name), h('td', { class: 'seq' }, p.seq), h('td', {}, `${p.seq.length} nt`)));
    sheet.append(h('h2', {}, 'Primers'), t);
  }
  sheet.append(h('h2', {}, 'Sequence'), h('pre', { class: 'seq' }, numbered(doc.seq)));
  return sheet;
}
