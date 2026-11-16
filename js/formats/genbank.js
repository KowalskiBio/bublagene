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
    if (!r) return [];
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

export function readGenBank(text, fileName = 'sequence.gb') {
  const records = text.split(/^\/\/\s*$/m).filter((r) => /LOCUS|ORIGIN/.test(r));
  if (!records.length) throw new Error('No GenBank record found');
  return records.map((rec) => readRecord(rec, fileName));
}

function readRecord(rec, fileName) {
  const lines = rec.split(/\r?\n/);
  const locus = lines.find((l) => l.startsWith('LOCUS')) ?? '';
  const name = locus.split(/\s+/)[1] || fileName.replace(/\.(gb|gbk|genbank|gbff)$/i, '');
  const circular = /\bcircular\b/i.test(locus);
  const doc = newDocument({ name, circular });
  doc.history = [{ when: new Date().toISOString(), action: 'imported', detail: `Imported from GenBank file ${fileName}` }];
  const def = [];
  let section = '';
  const featLines = [];
  let seq = '';
  for (const line of lines) {
    if (/^[A-Z]/.test(line)) section = line.split(/\s+/)[0];
    if (section === 'DEFINITION') def.push(line.replace(/^DEFINITION\s*/, '').trim());
    else if (section === 'SOURCE' && line.startsWith('  ORGANISM')) doc.notes.organism = line.replace(/^\s+ORGANISM\s+/, '').trim();
    else if (section === 'FEATURES' && !line.startsWith('FEATURES')) featLines.push(line);
    else if (section === 'ORIGIN' && !line.startsWith('ORIGIN')) seq += line;
  }
  doc.seq = cleanSequence(seq);
  doc.notes.description = def.join(' ').replace(/\.$/, '');
  const L = doc.seq.length;

  // group feature lines
  let cur = null;
  const feats = [];
  for (const line of featLines) {
    if (/^ {5}\S/.test(line)) {
      cur = { key: line.slice(5, 21).trim(), loc: line.slice(21).trim(), quals: [] };
      feats.push(cur);
    } else if (cur && /^ {21}/.test(line)) {
      const body = line.slice(21);
      if (body.startsWith('/')) cur.quals.push(body.slice(1));
      else if (cur.quals.length) cur.quals[cur.quals.length - 1] += (cur.quals[cur.quals.length - 1].endsWith('-') ? '' : ' ') + body.trim();
      else cur.loc += body.trim();
    }
  }
  for (const f of feats) {
    if (f.key === 'source') continue;
    const loc = parseLocationDetailed(f.loc, L);
    if (!loc) continue;
    const qualifiers = parseQualifiers(f.quals);
    const label = qualifiers.find((q) => ['label', 'gene', 'product', 'note'].includes(q.name));
    const color = colorFor(f.key);
    const feature = {
      id: newId(), name: label ? label.values[0].value.slice(0, 40) : f.key, type: f.key, strand: loc.strand,
      color, segments: loc.segments.map((s) => ({ ...s, color })),
      qualifiers, extra: {},
    };
    if (loc.locationType === 'order') feature.locationType = 'order';
    doc.features.push(feature);
  }
  doc.source = { format: 'genbank', fileName };
  doc.seqEdited = false;
  return doc;
}

// Qualifier lines -> [{name, values}]. Consecutive occurrences of one name are
// kept as one entry with several values; the order of first appearance is kept.
export function parseQualifiers(lines) {
  const out = [];
  for (const q of lines) {
    const eq = q.indexOf('=');
    const qn = eq < 0 ? q : q.slice(0, eq);
    let v = eq < 0 ? '' : q.slice(eq + 1);
    if (v.startsWith('"')) v = v.slice(1, v.endsWith('"') ? -1 : undefined).replace(/""/g, '"');
    if (qn === 'translation') v = v.replace(/\s+/g, '');
    const kind = /^-?\d+$/.test(v) ? 'int' : 'text';
    const prev = out[out.length - 1];
    if (prev && prev.name === qn) prev.values.push({ kind, value: v });
    else out.push({ name: qn, values: [{ kind, value: v }] });
  }
  return out;
}

function wrapQualifier(name, value) {
  const isNum = /^-?\d+$/.test(value);
  const text = isNum ? `/${name}=${value}` : `/${name}="${value.replace(/"/g, '""')}"`;
  const out = [];
  let rest = text;
  while (rest.length > 58) {
    let cut = rest.lastIndexOf(' ', 58);
    if (cut < 20) cut = 58;
    out.push(rest.slice(0, cut));
    rest = rest.slice(cut).trimStart();
  }
  out.push(rest);
  return out.map((l) => ' '.repeat(21) + l).join('\n');
}

export function writeGenBank(doc) {
  const L = doc.seq.length;
  const d = new Date();
  const mon = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'][d.getMonth()];
  const date = `${String(d.getDate()).padStart(2, '0')}-${mon}-${d.getFullYear()}`;
  const locusName = doc.name.replace(/\s+/g, '_').slice(0, 16).padEnd(16);
  const lines = [];
  lines.push(`LOCUS       ${locusName} ${String(L).padStart(11)} bp    DNA     ${doc.circular ? 'circular' : 'linear  '} SYN ${date}`);
  lines.push(`DEFINITION  ${doc.notes.description || doc.name}.`);
  lines.push('ACCESSION   .');
  lines.push('VERSION     .');
  lines.push('KEYWORDS    BublaGene.');
  lines.push(`SOURCE      ${doc.notes.organism || 'synthetic DNA construct'}`);
  lines.push(`  ORGANISM  ${doc.notes.organism || 'synthetic DNA construct'}`);
  lines.push('FEATURES             Location/Qualifiers');
  lines.push(`     source          1..${L}`);
  lines.push(wrapQualifier('mol_type', 'other DNA'));
  for (const f of doc.features) {
    lines.push(`     ${f.type.padEnd(16)}${formatLocation(f, L)}`);
    const quals = (f.qualifiers || []).map((q) => ({ name: q.name, values: [...q.values] }));
    if (f.name && f.name !== f.type && !quals.some((q) => ['label', 'gene', 'product', 'note'].includes(q.name))) {
      quals.unshift({ name: 'label', values: [{ value: f.name }] });
    }
    for (const q of quals) for (const v of q.values) lines.push(wrapQualifier(q.name, String(v.value ?? '')));
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

export function writeFasta(doc) {
  const lines = [`>${doc.name.replace(/\s+/g, '_')} ${doc.circular ? 'circular' : 'linear'} ${doc.seq.length} bp`];
  for (let i = 0; i < doc.seq.length; i += 70) lines.push(doc.seq.slice(i, i + 70));
  return lines.join('\n') + '\n';
}

export { rangeLength, translate, featureTableId, recomputeFeature };
