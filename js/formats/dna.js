// Binary .dna reader/writer.
//
// Reading keeps every original segment so that writing can pass through the
// parts BublaGene does not edit (history, display options, alignments, ...).

import { parseXML, child, children, el, toXML, htmlToText } from '../xml.js';
import { newDocument, newId, colorFor } from '../model.js';
import { mod } from '../seq.js';

export const SEGMENT_NAMES = [
  'Uncompressed DNA', 'Compressed DNA', 'Sequence Index', 'Restriction Sites Index',
  'Features (Binary)', 'Primers', 'Notes', 'History', 'Additional Sequence Properties',
  'File Description', 'Features (XML)', 'Ancestral Sequence', 'Additional Root Annotations',
  'Display Options', 'Local Enzyme Sets', 'Unsynced Case Changes', 'Alignable Sequence',
  'Alignable Sequences', 'Sequence Trace', 'Uracil Positions', 'Colors',
  'Uncompressed Protein', 'Numbering', 'Embedded Files', 'Alignment Summary', 'Alignment',
  'Edited Alignment', 'Compressed Alignment', 'Enzyme Visibility',
  'Ancestral Sequence Manipulation', 'Compressed Data', 'Data with ID', 'Uncompressed RNA',
  'Gel Summary', 'Secondary Structure', 'Crispr Sites',
];

const BASES = 'GATCNBDHKMRSVWY?';
const td = new TextDecoder('utf-8');
const te = new TextEncoder();
// Eight-byte header magic that identifies a binary .dna document.
const COOKIE_MAGIC = Uint8Array.from([0x53, 0x6e, 0x61, 0x70, 0x47, 0x65, 0x6e, 0x65]);

const u32 = (b, o) => ((b[o] << 24) | (b[o + 1] << 16) | (b[o + 2] << 8) | b[o + 3]) >>> 0;
const u16 = (b, o) => (b[o] << 8) | b[o + 1];

export function isXZ(d) {
  return d.length >= 6 && d[0] === 0xfd && d[1] === 0x37 && d[2] === 0x7a && d[3] === 0x58 && d[4] === 0x5a && d[5] === 0;
}

export function readSegments(bytes, start = 0, end = bytes.length) {
  const out = [];
  let i = start;
  while (i < end) {
    if (i + 5 > end) throw new Error(`Truncated segment header at byte ${i}`);
    const type = bytes[i];
    const n = u32(bytes, i + 1);
    if (i + 5 + n > end) throw new Error(`Segment 0x${type.toString(16)} at byte ${i} runs past the end of the file`);
    out.push({ type, data: bytes.subarray(i + 5, i + 5 + n) });
    i += 5 + n;
  }
  return out;
}

export function decodeCompressedDNA(d) {
  const seqLen = u32(d, 0);
  const flags = d[4];
  const blocks = u32(d, 5);
  let i = 13;
  const parts = [];
  for (let b = 0; b < blocks; b++) {
    const enc = d[i];
    const count = u32(d, i + 1);
    i += 5;
    const bits = enc === 1 ? 2 : enc === 2 ? 4 : 0;
    if (!bits) throw new Error(`Unknown compressed-DNA block encoding ${enc}`);
    const per = 8 / bits;
    const mask = (1 << bits) - 1;
    const nbytes = Math.ceil((count * bits) / 8);
    let s = '';
    for (let k = 0; k < nbytes && s.length < count; k++) {
      const byte = d[i + k];
      for (let j = 0; j < per && s.length < count; j++) s += BASES[(byte >> (8 - bits * (j + 1))) & mask];
    }
    parts.push(s);
    i += nbytes;
  }
  const seq = parts.join('');
  if (seq.length !== seqLen) throw new Error('Compressed DNA length mismatch');
  return { flags, seq };
}

// "1653-1655" (1-based inclusive) -> {start: 1652, end: 1655}; wraps when end < start
function parseRange(text, L) {
  const m = /^(\d+)-(\d+)$/.exec(text.trim());
  if (!m) {
    const n = parseInt(text, 10);
    return { start: n - 1, end: n };
  }
  const a = parseInt(m[1], 10);
  const b = parseInt(m[2], 10);
  return { start: a - 1, end: b <= L ? b : mod(b, L) || L };
}

function formatRange(r) {
  return `${r.start + 1}-${r.end}`;
}

function parseFeatures(xml, L) {
  const root = parseXML(xml);
  const out = [];
  for (const f of children(root, 'Feature')) {
    const { name, type, directionality, ...extra } = f.attrs;
    const strand = directionality === '1' ? 1 : directionality === '2' ? -1 : 0;
    const segments = children(f, 'Segment').map((s) => {
      const { range, color, name: segName, ...segExtra } = s.attrs;
      return { ...parseRange(range ?? '1-1', L), color, name: segName, extra: segExtra };
    });
    const qualifiers = children(f, 'Q').map((q) => ({
      name: q.attrs.name,
      values: children(q, 'V').map((v) => {
        if ('int' in v.attrs) return { kind: 'int', value: v.attrs.int };
        if ('predef' in v.attrs) return { kind: 'predef', value: v.attrs.predef };
        return { kind: 'text', value: htmlToText(v.attrs.text ?? v.text ?? '') };
      }),
    }));
    out.push({
      id: newId(), name: name ?? 'feature', type: type ?? 'misc_feature', strand,
      color: segments[0]?.color || colorFor(type), segments, qualifiers, extra,
    });
  }
  return { features: out, attrs: root.attrs };
}

function parsePrimers(xml) {
  const root = parseXML(xml);
  return {
    primers: children(root, 'Primer').map((p) => {
      const { name, sequence, description, ...extra } = p.attrs;
      return { id: newId(), name: name ?? 'primer', seq: (sequence ?? '').toUpperCase(), description: htmlToText(description ?? ''), extra };
    }),
    hybridization: child(root, 'HybridizationParams')?.attrs ?? null,
    attrs: root.attrs,
  };
}

const NOTE_FIELDS = {
  Description: 'description', Organism: 'organism', Comments: 'comments', Type: 'type',
  CreatedBy: 'createdBy', AccessionNumber: 'accession', SequenceClass: 'sequenceClass',
  TransformedInto: 'transformedInto', CustomMapLabel: 'customMapLabel', UUID: 'uuid',
};

function parseNotes(xml) {
  const root = parseXML(xml);
  const notes = {};
  for (const c of root.children) {
    if (NOTE_FIELDS[c.tag]) notes[NOTE_FIELDS[c.tag]] = htmlToText(c.text);
    if (c.tag === 'Created') notes.created = c.text;
    if (c.tag === 'LastModified') notes.modified = c.text;
    if (c.tag === 'References') {
      notes.references = children(c, 'Reference').map((r) => ({
        title: htmlToText(r.attrs.title ?? ''), journal: htmlToText(r.attrs.journal ?? ''),
        pubMedID: r.attrs.pubMedID, doi: r.attrs.doi,
      }));
    }
  }
  return { notes, root };
}

function parseProps(xml) {
  const root = parseXML(xml);
  const props = {};
  for (const c of root.children) props[c.tag] = c.text;
  return props;
}

export function readDna(buffer, fileName = 'sequence.dna') {
  const bytes = buffer instanceof Uint8Array ? buffer : new Uint8Array(buffer);
  const segs = readSegments(bytes);
  const cookie = segs.find((s) => s.type === 0x09);
  if (!cookie || !COOKIE_MAGIC.every((b, i) => cookie.data[i] === b)) {
    throw new Error('Not a supported .dna file (missing header)');
  }
  const docType = u16(cookie.data, 8);
  if (docType !== 1) throw new Error(`This .dna file is not a DNA document (type ${docType}); only DNA files are supported`);

  const doc = newDocument({ name: fileName.replace(/\.dna$/i, '') });
  doc.history = [];
  const warnings = [];
  const seqSeg = segs.find((s) => s.type === 0x00);
  const compSeg = segs.find((s) => s.type === 0x01);
  if (seqSeg) {
    doc.flags = seqSeg.data[0];
    doc.seq = td.decode(seqSeg.data.subarray(1)).toUpperCase();
  } else if (compSeg) {
    const r = decodeCompressedDNA(compSeg.data);
    doc.flags = r.flags;
    doc.seq = r.seq;
  } else {
    throw new Error('No DNA sequence segment found');
  }
  doc.circular = (doc.flags & 1) === 1;
  const L = doc.seq.length;

  const unparsed = [];
  const xmlOf = (type) => {
    const s = segs.find((x) => x.type === type);
    if (!s) return null;
    if (isXZ(s.data)) {
      warnings.push(`${SEGMENT_NAMES[type] ?? type} is xz-compressed; BublaGene can't read it and will save it unchanged`);
      unparsed.push(type);
      return null;
    }
    return td.decode(s.data);
  };

  const fx = xmlOf(0x0a);
  let featureAttrs = {};
  if (fx) ({ features: doc.features, attrs: featureAttrs } = parseFeatures(fx, L));
  const px = xmlOf(0x05);
  let primerInfo = { hybridization: null, attrs: {} };
  if (px) {
    primerInfo = parsePrimers(px);
    doc.primers = primerInfo.primers;
  }
  const nx = xmlOf(0x06);
  let notesRoot = null;
  if (nx) {
    const r = parseNotes(nx);
    doc.notes = { ...doc.notes, ...r.notes };
    notesRoot = r.root;
  }
  const ax = xmlOf(0x08);
  const props = ax ? parseProps(ax) : null;
  const historyCount = segs.filter((s) => s.type === 0x0b).length;
  if (historyCount) doc.history.push({ when: doc.notes.modified, action: 'imported', detail: `${historyCount} earlier versions stored in the file (kept, not shown)` });

  doc.source = {
    format: 'dna', fileName, cookie: cookie.data.slice(), segments: segs.map((s) => ({ type: s.type, data: s.data.slice() })),
    featureAttrs, primerInfo, notesRoot, props, warnings, unparsed,
  };
  doc.seqEdited = false;
  return doc;
}

// ---- writing ----------------------------------------------------------------

function segBytes(type, data) {
  const out = new Uint8Array(5 + data.length);
  out[0] = type;
  out[1] = (data.length >>> 24) & 255;
  out[2] = (data.length >>> 16) & 255;
  out[3] = (data.length >>> 8) & 255;
  out[4] = data.length & 255;
  out.set(data, 5);
  return out;
}

const html = (s) => (s ? `<html><body>${s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/\n/g, '<br>')}</body></html>` : s);

function featuresXML(doc) {
  const src = doc.source?.featureAttrs ?? {};
  const feats = doc.features.map((f, i) => {
    const attrs = { ...f.extra, recentID: f.extra?.recentID ?? 0, name: f.name, type: f.type };
    if (f.strand === 1) attrs.directionality = '1';
    else if (f.strand === -1) attrs.directionality = '2';
    const segsX = f.segments.map((s) => el('Segment', { ...s.extra, range: formatRange(s), color: s.color || f.color, name: s.name, type: s.extra?.type ?? 'standard' }));
    const qs = (f.qualifiers || []).map((q) => el('Q', { name: q.name }, q.values.map((v) => (v.kind === 'int' ? el('V', { int: v.value }) : v.kind === 'predef' ? el('V', { predef: v.value }) : el('V', { text: html(v.value) })))));
    return el('Feature', attrs, [...segsX, ...qs]);
  });
  return toXML(el('Features', { ...src, nextValidID: Math.max(Number(src.nextValidID) || 0, doc.features.length) }, feats), { declaration: true });
}

function primersXML(doc, bindingSitesFor) {
  const info = doc.source?.primerInfo ?? {};
  const hp = info.hybridization ?? { minContinuousMatchLen: '10', allowMismatch: '1', minMeltingTemperature: '40', showAdditionalFivePrimeMatches: '1', minimumFivePrimeAnnealing: '15' };
  const ps = doc.primers.map((p) => {
    const sites = bindingSitesFor ? bindingSitesFor(p) : [];
    const kids = sites.map((s) => {
      const loc = `${s.start + 1}-${s.end}`;
      return el('BindingSite', { location: loc, boundStrand: s.strand === 1 ? '0' : '1', annealedBases: s.annealed, meltingTemperature: Math.round(s.tm) },
        [el('Component', { hybridizedRange: loc, bases: s.annealed })]);
    });
    return el('Primer', { ...p.extra, recentID: p.extra?.recentID ?? 0, name: p.name, sequence: p.seq, description: html(p.description) }, kids);
  });
  return toXML(el('Primers', { ...(info.attrs ?? {}), nextValidID: Math.max(Number(info.attrs?.nextValidID) || 0, doc.primers.length) }, [el('HybridizationParams', hp), ...ps]), { declaration: true });
}

function snapDate(iso) {
  const d = iso ? new Date(iso) : new Date();
  if (Number.isNaN(d.getTime())) return { date: iso, utc: '0:0:0' };
  return { date: `${d.getUTCFullYear()}.${d.getUTCMonth() + 1}.${d.getUTCDate()}`, utc: `${d.getUTCHours()}:${d.getUTCMinutes()}:${d.getUTCSeconds()}` };
}

function uuid() {
  if (globalThis.crypto?.randomUUID) return globalThis.crypto.randomUUID();
  return 'xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx'.replace(/[xy]/g, (c) => ((Math.random() * 16) | (c === 'x' ? 0 : 8)).toString(16));
}

function notesXML(doc) {
  const n = doc.notes;
  const keepRoot = doc.source?.notesRoot;
  const set = new Map();
  const add = (tag, text, attrs = {}) => set.set(tag, el(tag, attrs, [], text ?? ''));
  add('UUID', n.uuid || uuid());
  add('Type', n.type || 'Synthetic');
  add('ConfirmedExperimentally', '0');
  if (!child(keepRoot, 'Created')) {
    const created = snapDate(n.created);
    add('Created', created.date, { UTC: created.utc });
  }
  const mod = snapDate(new Date().toISOString());
  add('LastModified', mod.date, { UTC: mod.utc });
  add('SequenceClass', n.sequenceClass || 'UNA');
  add('TransformedInto', n.transformedInto || 'unspecified');
  if (n.description) add('Description', html(n.description));
  if (n.organism) add('Organism', n.organism);
  if (n.comments) add('Comments', html(n.comments));
  if (n.createdBy) add('CreatedBy', n.createdBy);
  // keep anything else that was in the original Notes (references, map labels, ...)
  const kids = [];
  if (keepRoot) {
    for (const c of keepRoot.children) {
      if (set.has(c.tag)) { kids.push(set.get(c.tag)); set.delete(c.tag); } else kids.push(c);
    }
  }
  for (const v of set.values()) kids.push(v);
  return toXML(el('Notes', {}, kids));
}

function propsXML(doc) {
  const p = { UpstreamStickiness: '0', DownstreamStickiness: '0', UpstreamModification: 'FivePrimePhosphorylated', DownstreamModification: 'FivePrimePhosphorylated', ...(doc.source?.props ?? {}) };
  return toXML(el('AdditionalSequenceProperties', {}, Object.entries(p).map(([k, v]) => el(k, {}, [], v))));
}

// Segments regenerated from the model; everything else is passed through.
const REGENERATED = new Set([0x00, 0x01, 0x05, 0x06, 0x08, 0x09, 0x0a]);
// Caches the format rebuilds; dropped when the sequence changed so they can't go stale.
const SEQUENCE_CACHES = new Set([0x02, 0x03]);

export function writeDna(doc, { bindingSitesFor } = {}) {
  let flags = doc.flags ?? 0x1e;
  flags = doc.circular ? flags | 1 : flags & ~1;
  const cookie = doc.source?.cookie ?? Uint8Array.from([...COOKIE_MAGIC, 0, 1, 0, 15, 0, 19]);
  const parts = [segBytes(0x09, cookie)];
  const seqData = new Uint8Array(1 + doc.seq.length);
  seqData[0] = flags;
  seqData.set(te.encode(doc.seq), 1);
  parts.push(segBytes(0x00, seqData));
  const passThrough = new Set(doc.source?.unparsed ?? []);
  for (const s of doc.source?.segments ?? []) {
    if (REGENERATED.has(s.type) && !passThrough.has(s.type)) continue;
    if (s.type === 0x00 || s.type === 0x01 || s.type === 0x09) continue;
    if (doc.seqEdited && SEQUENCE_CACHES.has(s.type)) continue;
    parts.push(segBytes(s.type, s.data));
  }
  if (!passThrough.has(0x08)) parts.push(segBytes(0x08, te.encode(propsXML(doc))));
  if (!passThrough.has(0x0a)) parts.push(segBytes(0x0a, te.encode(featuresXML(doc))));
  if (!passThrough.has(0x05)) parts.push(segBytes(0x05, te.encode(primersXML(doc, bindingSitesFor))));
  if (!passThrough.has(0x06)) parts.push(segBytes(0x06, te.encode(notesXML(doc))));
  const total = parts.reduce((n, p) => n + p.length, 0);
  const out = new Uint8Array(total);
  let o = 0;
  for (const p of parts) { out.set(p, o); o += p.length; }
  return out;
}

export function describeSegments(doc) {
  return (doc.source?.segments ?? []).map((s) => ({ type: s.type, name: SEGMENT_NAMES[s.type] ?? `unknown 0x${s.type.toString(16)}`, bytes: s.data.length }));
}
