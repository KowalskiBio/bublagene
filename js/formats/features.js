// Feature-table formats: GFF3, BED, CSV and a plain-text map.

import { newDocument, newId, colorFor, featureSpan, qualifierText, setQualifier } from '../model.js';
import { cleanSequence } from '../seq.js';

const esc = (s) => String(s).replace(/([;=%&\t\n\r])/g, (c) => `%${c.charCodeAt(0).toString(16).toUpperCase().padStart(2, '0')}`);
const unesc = (s) => String(s).replace(/%([0-9A-Fa-f]{2})/g, (_, h) => String.fromCharCode(parseInt(h, 16)));
const attr = (f) => {
  const parts = [`ID=${esc(f.id)}`, `Name=${esc(f.name)}`];
  for (const q of f.qualifiers ?? []) {
    parts.push(`${esc(q.name)}=${esc(q.values.map((v) => v.value).join(','))}`);
  }
  return parts.join(';');
};

export function writeGFF3(doc) {
  const lines = ['##gff-version 3', `##sequence-region ${esc(doc.name)} 1 ${doc.seq.length}`];
  for (const f of doc.features ?? []) {
    const phase = f.type === 'CDS' ? Math.max(0, (Number(qualifierText(f, 'codon_start')) || 1) - 1) % 3 : '.';
    const strand = f.strand === -1 ? '-' : f.strand === 1 ? '+' : '.';
    for (const s of f.segments ?? []) {
      const start = Math.min(s.start, s.end) + 1;
      const end = Math.max(s.start, s.end);
      lines.push([esc(doc.name), 'BublaGene', esc(f.type), start, end, '.', strand, phase, attr(f)].join('\t'));
    }
  }
  return lines.join('\n') + '\n';
}

export function readGFF3(text, fileName = 'features.gff3') {
  const seqMatch = /##sequence-region\s+\S+\s+1\s+(\d+)/.exec(text);
  const length = seqMatch ? Number(seqMatch[1]) : 0;
  const doc = newDocument({ name: fileName.replace(/\.gff3?$/i, ''), seq: 'N'.repeat(length) });
  const groups = new Map();
  for (const line of text.split(/\r?\n/)) {
    if (!line || line.startsWith('#')) continue;
    const cols = line.split('\t');
    if (cols.length < 9) continue;
    const [seqid, source, type, start, end, score, strand, phase, attributes] = cols;
    const attrs = {};
    for (const pair of attributes.split(';')) {
      if (!pair) continue;
      const eq = pair.indexOf('=');
      if (eq < 0) continue;
      attrs[unesc(pair.slice(0, eq))] = unesc(pair.slice(eq + 1));
    }
    const id = attrs.ID ?? `gff-${groups.size + 1}-${start}`;
    if (!groups.has(id)) {
      groups.set(id, {
        id, name: attrs.Name ?? type, type, strand: strand === '-' ? -1 : strand === '+' ? 1 : 0,
        color: colorFor(type), segments: [], qualifiers: [], extra: {}, locationType: 'join',
      });
    }
    const f = groups.get(id);
    f.segments.push({ start: Number(start) - 1, end: Number(end), color: f.color });
    for (const [key, value] of Object.entries(attrs)) {
      if (key === 'ID' || key === 'Name') continue;
      f.qualifiers.push({ name: key, values: [{ kind: /^-?\d+$/.test(value) ? 'int' : 'text', value }] });
    }
    if (type === 'CDS' && phase !== '.' && phase !== '0') setQualifier(f, 'codon_start', String(Number(phase) + 1));
    f.extra.seqid = seqid;
  }
  for (const f of groups.values()) {
    f.segments.sort((a, b) => (f.strand === -1 ? b.start - a.start : a.start - b.start));
    doc.features.push(f);
  }
  doc.source = { format: 'gff3', fileName };
  return doc;
}

export function writeBED(doc) {
  const lines = ['track name="' + doc.name.replace(/"/g, '') + '"'];
  for (const f of doc.features ?? []) {
    for (const s of f.segments ?? []) {
      const start = Math.min(s.start, s.end);
      const end = Math.max(s.start, s.end);
      lines.push([doc.name, start, end, f.name.replace(/\s+/g, '_'), 0, f.strand === -1 ? '-' : f.strand === 1 ? '+' : '.'].join('\t'));
    }
  }
  return lines.join('\n') + '\n';
}

export function readBED(text, fileName = 'features.bed') {
  const doc = newDocument({ name: fileName.replace(/\.bed$/i, '') });
  let max = 0;
  const groups = new Map();
  for (const line of text.split(/\r?\n/)) {
    if (!line || line.startsWith('track') || line.startsWith('#')) continue;
    const cols = line.split('\t');
    if (cols.length < 3) continue;
    const start = Number(cols[1]);
    const end = Number(cols[2]);
    if (!Number.isFinite(start) || !Number.isFinite(end)) continue;
    max = Math.max(max, end);
    const name = cols[3] ?? 'feature';
    const strand = cols[5] === '-' ? -1 : cols[5] === '+' ? 1 : 0;
    const key = `${name}:${strand}`;
    if (!groups.has(key)) groups.set(key, { id: newId(), name, type: 'misc_feature', strand, color: colorFor('misc_feature'), segments: [], qualifiers: [], extra: {} });
    groups.get(key).segments.push({ start, end, color: colorFor('misc_feature') });
  }
  doc.seq = 'N'.repeat(max);
  doc.features = [...groups.values()];
  doc.source = { format: 'bed', fileName };
  return doc;
}

const csvCell = (v) => {
  const s = String(v ?? '');
  return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
};

export function writeFeatureCSV(doc) {
  const rows = ['name,type,strand,start,end,segments,qualifiers'];
  for (const f of doc.features ?? []) {
    const qualifiers = (f.qualifiers ?? []).map((q) => `${q.name}=${q.values.map((v) => v.value).join('|')}`).join('; ');
    const span = featureSpan(f);
    rows.push([f.name, f.type, f.strand, span.start + 1, span.end, f.segments.length, qualifiers].map(csvCell).join(','));
  }
  return rows.join('\n') + '\n';
}

function parseCSV(text) {
  const rows = [];
  let row = [];
  let cell = '';
  let quoted = false;
  for (let i = 0; i < text.length; i++) {
    const c = text[i];
    if (quoted) {
      if (c === '"' && text[i + 1] === '"') { cell += '"'; i++; }
      else if (c === '"') quoted = false;
      else cell += c;
    } else if (c === '"') quoted = true;
    else if (c === ',') { row.push(cell); cell = ''; }
    else if (c === '\n') { row.push(cell); rows.push(row); row = []; cell = ''; }
    else if (c !== '\r') cell += c;
  }
  if (cell || row.length) { row.push(cell); rows.push(row); }
  return rows.filter((r) => r.some((x) => x !== ''));
}

export function readFeatureCSV(text, fileName = 'features.csv') {
  const rows = parseCSV(text);
  if (!rows.length) throw new Error('No CSV rows found');
  const header = rows[0].map((h) => h.trim());
  const doc = newDocument({ name: fileName.replace(/\.csv$/i, '') });
  let max = 0;
  for (const row of rows.slice(1)) {
    const get = (name) => row[header.indexOf(name)] ?? '';
    const start = Number(get('start'));
    const end = Number(get('end'));
    if (!Number.isFinite(start) || !Number.isFinite(end)) continue;
    max = Math.max(max, end);
    const qualifiers = String(get('qualifiers')).split(';').filter(Boolean).map((part) => {
      const eq = part.indexOf('=');
      const name = eq < 0 ? part.trim() : part.slice(0, eq).trim();
      const value = eq < 0 ? '' : part.slice(eq + 1).trim();
      return { name, values: value.split('|').map((v) => ({ kind: 'text', value: v })) };
    });
    doc.features.push({
      id: newId(), name: get('name') || 'feature', type: get('type') || 'misc_feature',
      strand: Number(get('strand')) || 0, color: colorFor(get('type') || 'misc_feature'),
      segments: [{ start: start - 1, end, color: colorFor(get('type') || 'misc_feature') }],
      qualifiers, extra: {},
    });
  }
  doc.seq = 'N'.repeat(max);
  doc.source = { format: 'csv', fileName };
  return doc;
}

export function writePlainMap(doc, { width = 70 } = {}) {
  const L = doc.seq.length;
  const lines = [`${doc.name}  ${L} bp  ${doc.circular ? 'circular' : 'linear'}`];
  const step = Math.max(1, Math.ceil(L / width));
  const ruler = Array(width).fill(' ');
  for (let i = 0; i < L; i += step) {
    const p = Math.floor(i / step);
    const label = String(i + 1);
    for (let k = 0; k < label.length && p + k < width; k++) ruler[p + k] = label[k];
  }
  lines.push(' ' + ruler.join(''));
  for (const f of doc.features ?? []) {
    if (f.strand === -1) continue;
    const span = featureSpan(f);
    const row = Array(width).fill(' ');
    const a = Math.min(width - 1, Math.floor(span.start / step));
    const b = Math.min(width - 1, Math.max(a, Math.floor((span.end - 1) / step)));
    for (let p = a; p <= b; p++) row[p] = '=';
    row[a] = '>';
    lines.push(`${f.strand === 1 ? '>' : ' '} ${row.join('')} ${f.name}`);
  }
  for (const f of doc.features ?? []) {
    if (f.strand !== -1) continue;
    const span = featureSpan(f);
    const row = Array(width).fill(' ');
    const a = Math.min(width - 1, Math.floor(span.start / step));
    const b = Math.min(width - 1, Math.max(a, Math.floor((span.end - 1) / step)));
    for (let p = a; p <= b; p++) row[p] = '=';
    row[b] = '<';
    lines.push(`  ${row.join('')} ${f.name}`);
  }
  return lines.join('\n') + '\n';
}

// A FASTA/GenBank sequence is not implied by a feature table; callers that
// need bases should merge the table into their own document.
export function mergeFeatureTable(doc, table) {
  const L = doc.seq.length;
  for (const f of table.features ?? []) {
    const segs = f.segments.map((s) => ({ ...s, start: Math.min(s.start, L), end: Math.min(s.end, L) })).filter((s) => s.end > s.start);
    if (segs.length) doc.features.push({ ...f, segments: segs });
  }
  return doc;
}

export { cleanSequence };
