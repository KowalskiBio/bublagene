// GenBank flat-file reader/writer.
//
// Locations support plain, `<`/`>` fuzzy, `complement`, `join` and `order`, and
// wrapping joins across the origin. Qualifiers are kept in file order; repeated
// qualifier names become one entry with several values (never collapsed into a
// string). Unknown qualifier keys survive a round trip.

import { newDocument, newId, colorFor, featureTableId, recomputeFeature } from '../model.js';
import { cleanSequence, rangeLength, translate } from '../seq.js';

// Parse a GenBank location into {strand, segments} (0-based, end-exclusive).
export function parseLocation(loc, L) {
  const d = parseLocationDetailed(loc, L);
  return d ? { strand: d.strand, segments: d.segments } : null;
}

// Full form, including whether the top-level operator was join or order.
export function parseLocationDetailed(loc, L) {
  loc = String(loc).replace(/\s+/g, '');
  let strand = 1;
  const segs = [];
  let locationType = 'plain';

  const walk = (s, comp) => {
    if (s.startsWith('complement(') && s.endsWith(')')) return walk(s.slice(11, -1), !comp);
    const m = /^(join|order)\((.*)\)$/.exec(s);
    if (m) {
      if (m[1] === 'order') locationType = 'order';
      else if (locationType === 'plain') locationType = 'join';
      const parts = splitTop(m[2]);
      const out = parts.flatMap((p) => walk(p, false));
      return comp ? out.reverse().map((x) => ({ ...x, comp: !x.comp })) : out;
    }
    if (s.includes(':')) return []; // remote reference
    const between = /^(<?)(\d+)\^(\d+)$/.exec(s);
    if (between) {
      const a = parseInt(between[2], 10);
      return [{ start: a - 1, end: a, comp, partialStart: between[1] === '<', partialEnd: false, between: true }];
    }
    const r = /^(<?|>?)(\d+)(?:\.\.(<?|>?)(\d+))?$/.exec(s);
    if (!r) {
      // a position known only to lie within a range, "(3.5)..8": take the widest reading
      const w = /^([<>]?)(?:\((\d+)\.(\d+)\)|(\d+))(?:\.\.([<>]?)(?:\((\d+)\.(\d+)\)|(\d+)))?$/.exec(s);
      if (!w) return [];
      const from = parseInt(w[2] ?? w[4], 10);
      const hasEnd = w[6] !== undefined || w[8] !== undefined;
      const to = hasEnd ? parseInt(w[7] ?? w[8], 10) : w[3] !== undefined ? parseInt(w[3], 10) : from;
      return [{ start: from - 1, end: to, comp, partialStart: w[1] === '<', partialEnd: w[5] === '>' }];
    }
    const a = parseInt(r[2], 10);
    const partialStart = r[1] !== '';
    const b = r[4] ? parseInt(r[4], 10) : a;
    const partialEnd = r[3] === '>' || (r[1] === '>' && !r[4]) || /^<?>\d+$/.test(s) && !r[4];
    return [{ start: a - 1, end: b, comp, partialStart, partialEnd }];
  };

  const raw = walk(loc, false);
  if (!raw.length) return null;
  if (raw.every((x) => x.comp)) strand = -1;
  for (const x of raw) {
    segs.push({
      start: x.start,
      end: x.end <= L ? x.end : L,
      ...(x.partialStart ? { partialStart: true } : {}),
      ...(x.partialEnd ? { partialEnd: true } : {}),
      ...(x.between ? { between: true } : {}),
    });
  }
  // store in top-strand order
  if (strand === -1) segs.reverse();
  // merge a "join(N..L,1..M)" across the origin into one wrapping segment
  if (segs.length >= 2) {
    const last = segs[segs.length - 1];
    const first = segs[0];
    if (strand === 1 && last.start === 0 && segs[segs.length - 2]?.end === L && segs.length === 2) {
      return { strand, segments: [{ start: first.start, end: last.end }], locationType };
    }
    if (strand === -1 && first.end === L && segs[1]?.start === 0 && segs.length === 2) {
      return { strand, segments: [{ start: first.start, end: segs[1].end }], locationType };
    }
  }
  return { strand, segments: segs, locationType };
}

function splitTop(s) {
  const out = [];
  let depth = 0;
  let cur = '';
  for (const ch of s) {
    if (ch === '(') depth++;
    if (ch === ')') depth--;
    if (ch === ',' && depth === 0) { out.push(cur); cur = ''; } else cur += ch;
  }
  if (cur) out.push(cur);
  return out;
}

export function formatLocation(f, L) {
  const parts = [];
  for (const s of f.segments) {
    const left = s.partialStart ? '<' : '';
    const right = s.partialEnd ? '>' : '';
    if (s.between) parts.push(`${left}${s.start + 1}^${s.end + 1}`);
    else if (s.end > s.start) parts.push(s.end - s.start === 1 ? `${left}${s.start + 1}${right}` : `${left}${s.start + 1}..${right}${s.end}`);
    else { parts.push(`${left}${s.start + 1}${right}..${L}`); parts.push(`1..${s.end}`); }
  }
  const op = f.locationType === 'order' ? 'order' : 'join';
  let loc = parts.length > 1 ? `${op}(${parts.join(',')})` : parts[0] ?? '1..0';
  if (f.strand === -1) loc = `complement(${loc})`;
  return loc;
}

// ---- reading --------------------------------------------------------------------

// The qualifier that names a feature, in the order other tools use. A long /note
// is a description, not a name.
const NAME_QUALIFIERS = ['label', 'gene', 'product', 'standard_name', 'locus_tag', 'note'];

export function nameFromQualifiers(qualifiers, fallback) {
  for (const n of NAME_QUALIFIERS) {
    const v = qualifiers.find((q) => q.name === n)?.values[0]?.value;
    if (v && (n !== 'note' || String(v).length <= 60)) return String(v);
  }
  return fallback;
}

// Keys some programs write that are not in the INSDC feature table.
const KEY_ALIASES = { feature: 'misc_feature', primer: 'primer_bind' };

const HEX_COLOR = /^#[0-9a-f]{6}$/i;
const COLOR_QUALIFIERS = ['ApEinfo_fwdcolor', 'ApEinfo_revcolor'];

// The colour a file chose for a feature (ApE, SnapGene and others write these).
function colorFromQualifiers(qualifiers, strand) {
  const order = strand === -1 ? [...COLOR_QUALIFIERS].reverse() : COLOR_QUALIFIERS;
  for (const n of order) {
    const v = qualifiers.find((q) => q.name === n)?.values[0]?.value;
    if (HEX_COLOR.test(String(v))) return String(v);
  }
  return null;
}

export function readGenBank(text, fileName = 'sequence.gb') {
  // a byte-order mark or old line endings must not hide the LOCUS line
  const clean = String(text).replace(/^﻿/, '').replace(/\r\n?/g, '\n');
  const records = clean.split(/^\/\/\s*$/m).filter((r) => /^\s*(LOCUS|ORIGIN)/m.test(r));
  if (!records.length) throw new Error('No GenBank record found');
  return records.map((rec) => readRecord(rec, fileName));
}

// LOCUS name, length, molecule type, topology, division and date, in the layouts
// NCBI, SnapGene, ApE and others write.
function parseLocusLine(line) {
  const info = { name: '', moleculeType: '', circular: false, division: '' };
  const m = /^LOCUS\s+(\S+)\s+(\d+)\s+(?:bp|aa)\b\s*(.*)$/i.exec(line);
  if (!m) {
    info.name = line.split(/\s+/)[1] || '';
    info.circular = /\bcircular\b/i.test(line);
    return info;
  }
  info.name = m[1];
  const rest = m[3].trim().split(/\s+/).filter(Boolean);
  const at = rest.findIndex((t) => /^(linear|circular)$/i.test(t));
  info.circular = at >= 0 && /circular/i.test(rest[at]);
  info.moleculeType = (at >= 0 ? rest.slice(0, at) : rest.slice(0, 1)).join(' ');
  info.division = (at >= 0 ? rest.slice(at + 1) : []).find((t) => /^[A-Z]{3}$/.test(t)) ?? '';
  return info;
}

// Header sections the program does not edit (accession, references, comments...)
// are kept as written and put back when the file is saved.
const OWN_SECTIONS = new Set(['LOCUS', 'DEFINITION', 'FEATURES', 'ORIGIN', 'BASE', 'CONTIG']);

function readRecord(rec, fileName) {
  const lines = rec.split('\n');
  const locus = parseLocusLine(lines.find((l) => l.startsWith('LOCUS')) ?? '');
  const name = locus.name || fileName.replace(/\.(gb|gbk|genbank|gbff)$/i, '');
  const doc = newDocument({ name, circular: locus.circular });
  doc.history = [{ when: new Date().toISOString(), action: 'imported', detail: `Imported from GenBank file ${fileName}` }];
  const def = [];
  const blocks = [];
  let section = '';
  let block = null;
  const featLines = [];
  let seq = '';
  for (const line of lines) {
    if (/^[A-Z]/.test(line)) {
      section = line.split(/\s+/)[0];
      block = OWN_SECTIONS.has(section) ? null : { key: section, lines: [] };
      if (block) blocks.push(block);
    }
    if (block) block.lines.push(line);
    if (section === 'DEFINITION') def.push(line.replace(/^DEFINITION\s*/, '').trim());
    else if (section === 'SOURCE' && line.startsWith('  ORGANISM')) doc.notes.organism = line.replace(/^\s+ORGANISM\s+/, '').trim();
    else if (section === 'FEATURES' && !line.startsWith('FEATURES')) featLines.push(line);
    else if (section === 'ORIGIN' && !line.startsWith('ORIGIN')) seq += line;
  }
  doc.seq = cleanSequence(seq);
  const description = def.join(' ').replace(/\s+/g, ' ').trim().replace(/\.$/, '').trim();
  doc.notes.description = description;
  const L = doc.seq.length;

  // group feature lines; a "/" at the start of a line only begins a new qualifier
  // when the previous one is not still inside its quotes ("and/or" can wrap there)
  const open = (s) => (s.match(/"/g)?.length ?? 0) % 2 === 1;
  let cur = null;
  const feats = [];
  for (const line of featLines) {
    if (/^ {5}\S/.test(line)) {
      const key = line.slice(5, 21).trim();
      cur = { key, loc: line.slice(21).trim(), quals: [] };
      feats.push(cur);
    } else if (cur && /^ {21}/.test(line)) {
      const body = line.slice(21);
      const last = cur.quals[cur.quals.length - 1];
      if (body.startsWith('/') && !(last !== undefined && open(last))) cur.quals.push(body.slice(1));
      else if (last !== undefined) cur.quals[cur.quals.length - 1] = last + (last.endsWith('-') ? '' : ' ') + body.trim();
      else cur.loc += body.trim();
    }
  }
  const skipped = [];
  let sourceQualifiers = null;
  for (const f of feats) {
    if (f.key === 'source') { sourceQualifiers = parseQualifiers(f.quals); continue; }
    const loc = parseLocationDetailed(f.loc, L);
    if (!loc) { skipped.push(`${f.key} ${f.loc}`); continue; }
    const qualifiers = parseQualifiers(f.quals);
    const type = KEY_ALIASES[f.key] ?? f.key;
    const color = colorFromQualifiers(qualifiers, loc.strand) ?? colorFor(type);
    const feature = {
      id: newId(), name: nameFromQualifiers(qualifiers, type), type, strand: loc.strand,
      color, segments: loc.segments.map((s) => ({ ...s, color })),
      qualifiers, extra: {},
    };
    if (loc.locationType === 'order') feature.locationType = 'order';
    doc.features.push(feature);
  }
  doc.notes.genbank = {
    moleculeType: locus.moleculeType, division: locus.division, length: L, blocks: blocks.map((b) => b.lines),
    sourceQualifiers,
  };
  doc.source = { format: 'genbank', fileName, ...(skipped.length ? { warnings: skipped } : {}) };
  doc.seqEdited = false;
  return doc;
}

// Qualifier lines -> [{name, values}]. Consecutive occurrences of one name are
// kept as one entry with several values; the order of first appearance is kept.
// A qualifier with no "=" is a flag (/pseudo) and stays one.
export function parseQualifiers(lines) {
  const out = [];
  for (const q of lines) {
    const eq = q.indexOf('=');
    const qn = eq < 0 ? q.trim() : q.slice(0, eq);
    let v = eq < 0 ? '' : q.slice(eq + 1);
    const quoted = v.startsWith('"');
    if (quoted) v = v.slice(1, v.endsWith('"') ? -1 : undefined).replace(/""/g, '"');
    if (qn === 'translation') v = v.replace(/\s+/g, '');
    const kind = eq < 0 ? 'flag' : !quoted && /^-?\d+$/.test(v) ? 'int' : 'text';
    const prev = out[out.length - 1];
    if (prev && prev.name === qn) prev.values.push({ kind, value: v });
    else out.push({ name: qn, values: [{ kind, value: v }] });
  }
  return out;
}

// ---- writing --------------------------------------------------------------------

const LINE = 79;

// Wrap at spaces so a reader can put the value back together. A word longer than
// the line is left whole, past the margin, rather than cut and given a stray space.
// A translation has no spaces, so it is cut at the margin.
function wrapQualifier(name, v) {
  const value = String(v.value ?? '');
  const flag = v.kind === 'flag' || value === '';
  const numeric = v.kind === 'int' || (v.kind === undefined && /^-?\d+$/.test(value));
  const text = flag ? `/${name}` : numeric ? `/${name}=${value}` : `/${name}="${value.replace(/"/g, '""')}"`;
  const width = LINE - 21;
  const out = [];
  let rest = text;
  while (rest.length > width) {
    let cut;
    let skip = 1;
    if (name === 'translation') { cut = width; skip = 0; } else {
      cut = rest.lastIndexOf(' ', width);
      if (cut < 1) cut = rest.indexOf(' ', width);
      if (cut < 1) break;
    }
    out.push(rest.slice(0, cut));
    rest = rest.slice(cut + skip);
  }
  out.push(rest);
  return out.map((l) => ' '.repeat(21) + l);
}

function wrapWords(label, text) {
  const indent = 12; // header text starts in column 13
  const out = [];
  let line = '';
  for (const word of text.split(/\s+/).filter(Boolean)) {
    if (line && indent + line.length + 1 + word.length > LINE) { out.push(line); line = word; } else line = line ? `${line} ${word}` : word;
  }
  out.push(line);
  return out.map((l, i) => (i === 0 ? label.padEnd(indent) : ' '.repeat(indent)) + l);
}

// Locations that do not fit on one line continue after a comma, as other tools write them.
function featureLines(type, loc) {
  const head = `${type} `.padEnd(16);
  const first = `     ${head}`;
  if (first.length + loc.length <= LINE) return [first + loc];
  const out = [];
  let line = first;
  for (const piece of loc.split(/(?<=,)/)) {
    if (line.length + piece.length > LINE && line.trim() !== type) { out.push(line); line = ' '.repeat(21); }
    line += piece;
  }
  out.push(line);
  return out;
}

function setQualifier(quals, name, value) {
  const q = quals.find((x) => x.name === name);
  if (q) q.values[0] = { kind: 'text', value };
  else quals.push({ name, values: [{ kind: 'text', value }] });
}

export function writeGenBank(doc) {
  const L = doc.seq.length;
  const gb = doc.notes?.genbank ?? {};
  const d = new Date();
  const mon = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'][d.getMonth()];
  const date = `${String(d.getDate()).padStart(2, '0')}-${mon}-${d.getFullYear()}`;
  const mol = gb.moleculeType || 'DNA';
  const locusName = String(doc.name).trim().replace(/\s+/g, '_') || 'Untitled';
  const lines = [];
  lines.push(`LOCUS       ${locusName.padEnd(16)} ${String(L).padStart(11)} bp    ${`${mol} `.padEnd(8)}${(doc.circular ? 'circular' : 'linear').padEnd(8)} ${gb.division || 'SYN'} ${date}`);
  lines.push(...wrapWords('DEFINITION', doc.notes?.description ? `${doc.notes.description}.` : '.'));
  if (gb.blocks?.length) {
    // the rest of the header as it was read; a reference to "bases 1 to N" follows the sequence length
    for (const block of gb.blocks) {
      for (const l of block) lines.push(l.replace(/\(bases (\d+) to (\d+)\)/g, (m, a, b) => (Number(b) === gb.length ? `(bases ${a} to ${L})` : m)));
    }
  } else {
    lines.push('ACCESSION   .');
    lines.push('VERSION     .');
    lines.push('KEYWORDS    .');
    lines.push(`SOURCE      ${doc.notes?.organism || 'synthetic DNA construct'}`);
    lines.push(`  ORGANISM  ${doc.notes?.organism || 'synthetic DNA construct'}`);
  }
  lines.push('FEATURES             Location/Qualifiers');
  lines.push(`     source          1..${L}`);
  const sourceQuals = gb.sourceQualifiers?.length ? gb.sourceQualifiers : [{ name: 'mol_type', values: [{ kind: 'text', value: /RNA/i.test(mol) ? 'other RNA' : 'other DNA' }] }];
  for (const q of sourceQuals) for (const v of q.values) lines.push(...wrapQualifier(q.name, v));
  for (const f of doc.features) {
    lines.push(...featureLines(f.type, formatLocation(f, L)));
    const quals = (f.qualifiers || []).map((q) => ({ name: q.name, values: q.values.map((v) => ({ ...v })) }));
    // renaming a feature in the program must reach the file
    if (f.name && f.name !== nameFromQualifiers(quals, f.type)) {
      if (quals.some((q) => q.name === 'label')) setQualifier(quals, 'label', f.name);
      else quals.unshift({ name: 'label', values: [{ kind: 'text', value: f.name }] });
    }
    // so must its colour; keep the file clean when it is just the default for its type
    const hasColor = quals.some((q) => COLOR_QUALIFIERS.includes(q.name));
    if (HEX_COLOR.test(f.color ?? '') && (hasColor || f.color.toLowerCase() !== colorFor(f.type).toLowerCase())) {
      for (const n of COLOR_QUALIFIERS) setQualifier(quals, n, f.color);
    }
    for (const q of quals) for (const v of q.values) lines.push(...wrapQualifier(q.name, v));
  }
  lines.push('ORIGIN');
  for (let i = 0; i < L; i += 60) {
    const chunk = doc.seq.slice(i, i + 60).toLowerCase();
    const groups = chunk.match(/.{1,10}/g).join(' ');
    lines.push(`${String(i + 1).padStart(9)} ${groups}`);
  }
  lines.push('//');
  return lines.join('\n') + '\n';
}

export function readFasta(text, fileName = 'sequence.fasta') {
  const docs = [];
  let cur = null;
  for (const line of text.split(/\r?\n/)) {
    if (line.startsWith('>')) {
      cur = { name: line.slice(1).trim().split(/\s+/)[0] || 'sequence', desc: line.slice(1).trim(), seq: '' };
      docs.push(cur);
    } else if (line.trim()) {
      if (!cur) { cur = { name: fileName.replace(/\.\w+$/, ''), desc: '', seq: '' }; docs.push(cur); }
      cur.seq += line;
    }
  }
  if (!docs.length) throw new Error('No sequence found');
  return docs.map((r) => {
    const doc = newDocument({ name: r.name, seq: cleanSequence(r.seq), circular: /circular/i.test(r.desc) });
    doc.notes.description = r.desc;
    doc.history = [{ when: new Date().toISOString(), action: 'imported', detail: `Imported from FASTA file ${fileName}` }];
    doc.source = { format: 'fasta', fileName };
    doc.seqEdited = false;
    return doc;
  });
}

// Streaming FASTA: push chunks, get one document at a time. The record length
// cap prevents a malformed single record from using unbounded memory.
export function createFastaParser({ maxRecordLength = 1e9, fileName = 'sequence.fasta' } = {}) {
  let buffer = '';
  let header = null;
  let seq = '';
  let records = 0;
  const make = () => {
    if (header == null) return null;
    const doc = newDocument({ name: header.split(/\s+/)[0] || 'sequence', seq: cleanSequence(seq), circular: /circular/i.test(header) });
    doc.notes.description = header;
    doc.source = { format: 'fasta', fileName };
    doc.seqEdited = false;
    records++;
    return doc;
  };
  return {
    push(chunk) {
      buffer += typeof chunk === 'string' ? chunk : new TextDecoder().decode(chunk, { stream: true });
      const out = [];
      let nl;
      while ((nl = buffer.indexOf('\n')) >= 0) {
        let line = buffer.slice(0, nl);
        buffer = buffer.slice(nl + 1);
        if (line.endsWith('\r')) line = line.slice(0, -1);
        if (line.startsWith('>')) {
          const doc = make();
          if (doc) out.push(doc);
          header = line.slice(1).trim();
          seq = '';
        } else if (line.trim()) {
          if (header == null) header = fileName.replace(/\.\w+$/, '');
          seq += line.trim();
          if (seq.length > maxRecordLength) throw new Error(`FASTA record exceeds ${maxRecordLength} bases`);
        }
      }
      return out;
    },
    end() {
      const out = [];
      if (buffer) {
        let line = buffer.replace(/\r$/, '');
        if (line.startsWith('>')) {
          const doc = make();
          if (doc) out.push(doc);
          header = line.slice(1).trim();
          seq = '';
        } else if (line.trim()) {
          if (header == null) header = fileName.replace(/\.\w+$/, '');
          seq += line.trim();
        }
        buffer = '';
      }
      const doc = make();
      if (doc) out.push(doc);
      header = null;
      return out;
    },
    get count() { return records; },
  };
}

export async function* streamFasta(iterable, options = {}) {
  const parser = createFastaParser(options);
  for await (const chunk of iterable) {
    for (const doc of parser.push(chunk)) yield doc;
  }
  for (const doc of parser.end()) yield doc;
}

export function writeFasta(doc) {
  const lines = [`>${doc.name.replace(/\s+/g, '_')} ${doc.circular ? 'circular' : 'linear'} ${doc.seq.length} bp`];
  for (let i = 0; i < doc.seq.length; i += 70) lines.push(doc.seq.slice(i, i + 70));
  return lines.join('\n') + '\n';
}

export { rangeLength, translate, featureTableId, recomputeFeature };
