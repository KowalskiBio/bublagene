// Diagram exports: SVG, PNG and a small vector PDF.
//
// The SVG/PDF writers use vector primitives. PNG is rasterized here (no canvas
// dependency) and encoded as an uncompressed RGBA PNG, so exports are
// deterministic and work in Node as well as the browser.

import { featureSpan, colorFor } from './model.js';
import { revcomp } from './seq.js';

const hex = (c) => {
  const m = /^#?([0-9a-f]{6})$/i.exec(String(c ?? ''));
  if (!m) return [120, 140, 160];
  const n = parseInt(m[1], 16);
  return [(n >> 16) & 255, (n >> 8) & 255, n & 255];
};

function lanes(features) {
  const rows = [];
  for (const f of features) {
    const span = featureSpan(f);
    let lane = rows.findIndex((end) => span.start >= end + 4);
    if (lane < 0) { lane = rows.length; rows.push(0); }
    rows[lane] = span.end;
    f.__lane = lane;
  }
  return rows.length;
}

export function exportSVG(doc, { width = 900, circular = doc.circular } = {}) {
  return circular ? circularSVG(doc, width) : linearSVG(doc, width);
}

function linearSVG(doc, width) {
  const feats = doc.features ?? [];
  const nLanes = lanes(feats);
  const height = 60 + nLanes * 22 + 30;
  const pad = 30;
  const scale = (width - pad * 2) / Math.max(1, doc.seq.length);
  const x = (p) => pad + p * scale;
  const parts = [`<svg xmlns="http://www.w3.org/2000/svg" width="${width}" height="${height}" viewBox="0 0 ${width} ${height}">`];
  parts.push(`<rect width="100%" height="100%" fill="#ffffff"/>`);
  parts.push(`<text x="${pad}" y="22" font-family="sans-serif" font-size="14" fill="#222">${escapeXml(doc.name)} — ${doc.seq.length} bp</text>`);
  parts.push(`<line x1="${pad}" y1="40" x2="${width - pad}" y2="40" stroke="#667" stroke-width="1"/>`);
  const tickStep = Math.max(1, Math.pow(10, Math.floor(Math.log10(Math.max(1, doc.seq.length / 8)))));
  for (let p = 0; p <= doc.seq.length; p += tickStep) {
    parts.push(`<line x1="${x(p).toFixed(1)}" y1="36" x2="${x(p).toFixed(1)}" y2="44" stroke="#667"/>`);
    parts.push(`<text x="${x(p).toFixed(1)}" y="54" font-family="sans-serif" font-size="9" fill="#556" text-anchor="middle">${p + 1}</text>`);
  }
  for (const f of feats) {
    const span = featureSpan(f);
    const y = 64 + (f.__lane ?? 0) * 22;
    const x1 = x(span.start);
    const x2 = Math.max(x1 + 2, x(span.end));
    const color = f.color ?? colorFor(f.type);
    parts.push(`<rect x="${x1.toFixed(1)}" y="${y}" width="${(x2 - x1).toFixed(1)}" height="14" rx="2" fill="${color}" opacity="0.85"/>`);
    if (f.strand === -1) parts.push(`<path d="M ${x1.toFixed(1)} ${y + 7} l 7 -5 v 10 z" fill="#333"/>`);
    else if (f.strand === 1) parts.push(`<path d="M ${x2.toFixed(1)} ${y + 7} l -7 -5 v 10 z" fill="#333"/>`);
    parts.push(`<text x="${x1.toFixed(1)}" y="${y + 28}" font-family="sans-serif" font-size="9" fill="#333">${escapeXml(f.name)}</text>`);
  }
  parts.push('</svg>');
  return parts.join('\n');
}

function circularSVG(doc, width) {
  const L = doc.seq.length;
  const size = width;
  const cx = size / 2;
  const cy = size / 2;
  const radius = size * 0.36;
  const parts = [`<svg xmlns="http://www.w3.org/2000/svg" width="${size}" height="${size}" viewBox="0 0 ${size} ${size}">`];
  parts.push(`<rect width="100%" height="100%" fill="#ffffff"/>`);
  parts.push(`<circle cx="${cx}" cy="${cy}" r="${radius}" fill="none" stroke="#667" stroke-width="1.2"/>`);
  parts.push(`<text x="${cx}" y="22" text-anchor="middle" font-family="sans-serif" font-size="14" fill="#222">${escapeXml(doc.name)} — ${L} bp circular</text>`);
  for (const f of doc.features ?? []) {
    for (const s of f.segments ?? []) {
      let len = s.end > s.start ? s.end - s.start : L - s.start + s.end;
      if (len <= 0 || len >= L) len = L - 1;
      const a0 = (s.start / L) * Math.PI * 2 - Math.PI / 2;
      const a1 = ((s.start + len) / L) * Math.PI * 2 - Math.PI / 2;
      const p0 = [cx + radius * Math.cos(a0), cy + radius * Math.sin(a0)];
      const p1 = [cx + radius * Math.cos(a1), cy + radius * Math.sin(a1)];
      const large = a1 - a0 > Math.PI ? 1 : 0;
      parts.push(`<path d="M ${p0[0].toFixed(1)} ${p0[1].toFixed(1)} A ${radius} ${radius} 0 ${large} 1 ${p1[0].toFixed(1)} ${p1[1].toFixed(1)}" fill="none" stroke="${f.color ?? colorFor(f.type)}" stroke-width="8" opacity="0.8"/>`);
    }
  }
  const tick = Math.max(1, Math.pow(10, Math.floor(Math.log10(Math.max(1, L / 8)))));
  for (let p = 0; p < L; p += tick) {
    const a = (p / L) * Math.PI * 2 - Math.PI / 2;
    const r0 = radius + 8;
    const r1 = radius + 14;
    const t0 = [cx + r0 * Math.cos(a), cy + r0 * Math.sin(a)];
    const t1 = [cx + r1 * Math.cos(a), cy + r1 * Math.sin(a)];
    parts.push(`<line x1="${t0[0].toFixed(1)}" y1="${t0[1].toFixed(1)}" x2="${t1[0].toFixed(1)}" y2="${t1[1].toFixed(1)}" stroke="#667"/>`);
    parts.push(`<text x="${(cx + (r1 + 10) * Math.cos(a)).toFixed(1)}" y="${(cy + (r1 + 10) * Math.sin(a)).toFixed(1)}" font-family="sans-serif" font-size="9" fill="#556" text-anchor="middle">${p + 1}</text>`);
  }
  parts.push('</svg>');
  return parts.join('\n');
}

function escapeXml(s) {
  return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
}

// ---- PNG --------------------------------------------------------------------

const CRC_TABLE = (() => {
  const table = new Uint32Array(256);
  for (let n = 0; n < 256; n++) {
    let c = n;
    for (let k = 0; k < 8; k++) c = c & 1 ? 0xedb88320 ^ (c >>> 1) : c >>> 1;
    table[n] = c >>> 0;
  }
  return table;
})();

function crc32(bytes) {
  let c = 0xffffffff;
  for (const b of bytes) c = CRC_TABLE[(c ^ b) & 0xff] ^ (c >>> 8);
  return (c ^ 0xffffffff) >>> 0;
}

function adler32(bytes) {
  let a = 1;
  let b = 0;
  for (const x of bytes) { a = (a + x) % 65521; b = (b + a) % 65521; }
  return ((b << 16) | a) >>> 0;
}

function chunk(type, data) {
  const out = new Uint8Array(12 + data.length);
  const dv = new DataView(out.buffer);
  dv.setUint32(0, data.length, false);
  for (let i = 0; i < 4; i++) out[4 + i] = type.charCodeAt(i);
  out.set(data, 8);
  dv.setUint32(8 + data.length, crc32(out.subarray(4, 8 + data.length)), false);
  return out;
}

function zlibStore(raw) {
  const blocks = Math.max(1, Math.ceil(raw.length / 65535));
  const out = [0x78, 0x01];
  for (let b = 0; b < blocks; b++) {
    const slice = raw.subarray(b * 65535, Math.min(raw.length, (b + 1) * 65535));
    out.push(b === blocks - 1 ? 1 : 0);
    out.push(slice.length & 255, (slice.length >>> 8) & 255);
    out.push(~slice.length & 255, (~slice.length >>> 8) & 255);
    out.push(...slice);
  }
  const ad = adler32(raw);
  out.push((ad >>> 24) & 255, (ad >>> 16) & 255, (ad >>> 8) & 255, ad & 255);
  return Uint8Array.from(out);
}

export function exportPNG(doc, { width = 900, height = null } = {}) {
  const L = Math.max(1, doc.seq.length);
  const nLanes = lanes(doc.features ?? []);
  const h = height ?? Math.max(120, 50 + nLanes * 16);
  const bytes = new Uint8Array(width * h * 4).fill(255);
  const put = (x, y, [r, g, b], alpha = 255) => {
    if (x < 0 || y < 0 || x >= width || y >= h) return;
    const i = (y * width + x) * 4;
    bytes[i] = Math.round(bytes[i] * (1 - alpha / 255) + r * (alpha / 255));
    bytes[i + 1] = Math.round(bytes[i + 1] * (1 - alpha / 255) + g * (alpha / 255));
    bytes[i + 2] = Math.round(bytes[i + 2] * (1 - alpha / 255) + b * (alpha / 255));
  };
  const rect = (x1, y1, x2, y2, color) => {
    for (let y = y1; y < y2; y++) for (let x = x1; x < x2; x++) put(x, y, color);
  };
  rect(0, 0, width, h, [255, 255, 255]);
  const pad = 20;
  const scale = (width - pad * 2) / L;
  rect(pad, 30, width - pad, 32, [100, 110, 130]);
  for (const f of doc.features ?? []) {
    const span = featureSpan(f);
    const x1 = Math.round(pad + span.start * scale);
    const x2 = Math.max(x1 + 2, Math.round(pad + span.end * scale));
    const y = 40 + (f.__lane ?? 0) * 16;
    const color = hex(f.color ?? colorFor(f.type));
    rect(x1, y, x2, y + 10, color);
    const tip = f.strand === -1 ? x1 : x2;
    for (let k = 0; k < 6; k++) for (let yy = 0; yy < 10 - k; yy++) put(tip + (f.strand === -1 ? k : -k), y + yy, [40, 40, 40]);
  }
  const stride = width * 4 + 1;
  const raw = new Uint8Array(stride * h);
  for (let y = 0; y < h; y++) {
    raw[y * stride] = 0;
    raw.set(bytes.subarray(y * width * 4, (y + 1) * width * 4), y * stride + 1);
  }
  const sig = Uint8Array.from([137, 80, 78, 71, 13, 10, 26, 10]);
  const ihdr = new Uint8Array(13);
  const dv = new DataView(ihdr.buffer);
  dv.setUint32(0, width, false);
  dv.setUint32(4, h, false);
  ihdr[8] = 8; ihdr[9] = 6; ihdr[10] = 0; ihdr[11] = 0; ihdr[12] = 0;
  const chunks = [sig, chunk('IHDR', ihdr), chunk('IDAT', zlibStore(raw)), chunk('IEND', new Uint8Array(0))];
  const total = chunks.reduce((n, c) => n + c.length, 0);
  const out = new Uint8Array(total);
  let o = 0;
  for (const c of chunks) { out.set(c, o); o += c.length; }
  return out;
}

// ---- PDF --------------------------------------------------------------------

const pdfEscape = (s) => String(s).replace(/\\/g, '\\\\').replace(/\(/g, '\\(').replace(/\)/g, '\\)');

export function exportPDF(doc, { width = 842, height = 595 } = {}) {
  const feats = doc.features ?? [];
  const nLanes = lanes(feats);
  const pad = 40;
  const scale = (width - pad * 2) / Math.max(1, doc.seq.length);
  const ops = [];
  ops.push(`0.4 0.45 0.55 RG 0.8 w`);
  ops.push(`${pad} ${height - 60} m ${width - pad} ${height - 60} l S`);
  ops.push(`0 0 0 rg BT /F1 12 Tf ${pad} ${height - 30} Td (${pdfEscape(`${doc.name} - ${doc.seq.length} bp`)}) Tj ET`);
  for (const f of feats) {
    const span = featureSpan(f);
    const x = pad + span.start * scale;
    const w = Math.max(1, (span.end - span.start) * scale);
    const y = height - 72 - (f.__lane ?? 0) * 16;
    const [r, g, b] = hex(f.color ?? colorFor(f.type));
    ops.push(`${(r / 255).toFixed(3)} ${(g / 255).toFixed(3)} ${(b / 255).toFixed(3)} rg ${x.toFixed(2)} ${y.toFixed(2)} ${w.toFixed(2)} 10 re f`);
    ops.push(`0 0 0 rg BT /F1 7 Tf ${x.toFixed(2)} ${(y - 8).toFixed(2)} Td (${pdfEscape(f.name)}) Tj ET`);
  }
  const content = ops.join('\n');
  const objects = [
    '<< /Type /Catalog /Pages 2 0 R >>',
    '<< /Type /Pages /Kids [3 0 R] /Count 1 >>',
    `<< /Type /Page /Parent 2 0 R /MediaBox [0 0 ${width} ${height}] /Resources << /Font << /F1 5 0 R >> >> /Contents 4 0 R >>`,
    `<< /Length ${content.length} >>\nstream\n${content}\nendstream`,
    '<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>',
  ];
  let out = '%PDF-1.4\n';
  const offsets = [0];
  for (const [i, obj] of objects.entries()) {
    offsets.push(out.length);
    out += `${i + 1} 0 obj\n${obj}\nendobj\n`;
  }
  const xref = out.length;
  out += `xref\n0 ${objects.length + 1}\n0000000000 65535 f \n`;
  for (let i = 1; i <= objects.length; i++) out += `${String(offsets[i]).padStart(10, '0')} 00000 n \n`;
  out += `trailer\n<< /Size ${objects.length + 1} /Root 1 0 R >>\nstartxref\n${xref}\n%%EOF\n`;
  return out;
}

export function exportDiagram(doc, format, options = {}) {
  const f = String(format).toLowerCase();
  if (f === 'svg') return exportSVG(doc, options);
  if (f === 'png') return exportPNG(doc, options);
  if (f === 'pdf') return exportPDF(doc, options);
  throw new Error(`Unknown diagram format ${format}`);
}

export { revcomp };
