#!/usr/bin/env node
// BublaGene command line: a way to use the library before there is a UI.
//
//   bg info      <file>                  summary, features, primers
//   bg seq       <file> [start..end]     sequence (1-based, inclusive)
//   bg cut       <file> [enzyme,...]     restriction sites and fragment sizes
//   bg orfs      <file> [minCodons]      open reading frames
//   bg primers   <file> <start..end>     design a primer pair for a region
//   bg pcr       <file> <fwd> <rev>      simulate PCR, print the product
//   bg clone     <vector> <insert> <enzyme,...> [out]   restriction cloning
//   bg gibson    <f1> <f2> [...] [-o out] [--linear]    Gibson assembly
//   bg convert   <file> <out.gb|.fasta|.dna>
//
// Files may be .dna (binary), .gb/.gbk/.genbank, or .fasta/.fa.

import { readFileSync, writeFileSync } from 'node:fs';
import { basename, extname } from 'node:path';
import { readDna, writeDna, describeSegments } from '../js/formats/dna.js';
import { readGenBank, writeGenBank, readFasta, writeFasta } from '../js/formats/genbank.js';
import { ENZYMES, ENZYME_BY_NAME, COMMON_SET, findSites, digest, cutEnd, overhangLabel } from '../js/enzymes.js';
import { findBindingSites, designPrimers, simulatePCR } from '../js/primers.js';
import { restrictionClone, gibsonAssemble } from '../js/cloning.js';
import { codingSequence, featureLength, featureSpan, isCodingSegment } from '../js/model.js';
import { gcContent, meltingTemp, translate, findOrfs, subseq, rangeLength } from '../js/seq.js';

const [, , cmd, file, ...rest] = process.argv;

function load(path) {
  const name = basename(path);
  const ext = extname(path).toLowerCase();
  if (ext === '.dna') return readDna(readFileSync(path), name);
  const text = readFileSync(path, 'utf8');
  if (['.gb', '.gbk', '.genbank', '.gbff'].includes(ext)) return readGenBank(text, name)[0];
  if (['.fasta', '.fa', '.fna', '.fas'].includes(ext)) return readFasta(text, name)[0];
  if (text.startsWith('>')) return readFasta(text, name)[0];
  if (/^LOCUS/m.test(text)) return readGenBank(text, name)[0];
  throw new Error(`Don't know how to read ${name}`);
}

function parseRegion(spec, L) {
  const m = /^(\d+)\.\.(\d+)$/.exec(spec ?? '');
  if (!m) throw new Error('Give a region as start..end, 1-based and inclusive, e.g. 100..600');
  return { start: Number(m[1]) - 1, end: Number(m[2]) };
}

const pad = (s, n) => String(s).padEnd(n);
const num = (n) => Number(n).toLocaleString('en-US');

function info(doc) {
  console.log(`${doc.name}  ${num(doc.seq.length)} bp  ${doc.circular ? 'circular' : 'linear'}  GC ${gcContent(doc.seq).toFixed(1)}%`);
  if (doc.notes.description) console.log(doc.notes.description);
  if (doc.source?.warnings?.length) for (const w of doc.source.warnings) console.log(`note: ${w}`);
  const L = doc.seq.length;
  if (doc.features.length) {
    console.log(`\nFeatures (${doc.features.length})`);
    for (const f of doc.features) {
      const sp = featureSpan(f);
      const loc = `${sp.start + 1}..${sp.end}${f.segments.length > 1 ? ` (${f.segments.length} segments)` : ''}`;
      const strand = f.strand === 1 ? '+' : f.strand === -1 ? '-' : ' ';
      let extra = '';
      if (f.type === 'CDS') {
        const cds = codingSequence(doc, f);
        const prot = translate(cds);
        const flags = [];
        if (cds.length % 3) flags.push('not a whole number of codons');
        if (prot.slice(0, -1).includes('*')) flags.push('internal stop');
        extra = `  ${prot.length - (prot.endsWith('*') ? 1 : 0)} aa${flags.length ? ` (${flags.join(', ')})` : ''}`;
      }
      const introns = f.segments.filter((s) => !isCodingSegment(s)).length;
      const spanNote = introns ? `, ${introns} gap${introns > 1 ? 's' : ''}` : '';
      console.log(`  ${pad(f.type, 14)} ${strand} ${pad(loc.replace(/\)$/, spanNote + ')'), 26)} ${pad(num(featureLength(f, L)) + ' bp', 10)} ${f.name}${extra}`);
    }
  }
  if (doc.primers.length) {
    console.log(`\nPrimers (${doc.primers.length})`);
    for (const p of doc.primers) {
      const sites = findBindingSites(doc.seq, doc.circular, p.seq);
      const where = sites.map((s) => `${s.start + 1}..${s.end}${s.strand === 1 ? '+' : '-'}`).join(', ') || 'no binding site';
      console.log(`  ${pad(p.name, 18)} ${pad(p.seq.length + ' nt', 7)} Tm ${meltingTemp(p.seq).toFixed(1)}°C  ${where}`);
    }
  }
  const segs = describeSegments(doc);
  if (segs.length) {
    console.log('\nContainer segments in the file');
    for (const s of segs) console.log(`  0x${s.type.toString(16).padStart(2, '0')}  ${pad(s.name, 32)} ${num(s.bytes)} bytes`);
  }
}

function cut(doc, names) {
  const list = (names ? names.split(',') : COMMON_SET).map((n) => n.trim()).filter(Boolean);
  const unknown = list.filter((n) => !ENZYME_BY_NAME[n]);
  if (unknown.length) throw new Error(`Unknown enzyme(s): ${unknown.join(', ')}. ${ENZYMES.length} are available; try: ${COMMON_SET.slice(0, 8).join(', ')}…`);
  const rows = [];
  for (const n of list) {
    const sites = findSites(doc.seq, doc.circular, ENZYME_BY_NAME[n]);
    if (sites.length) rows.push({ n, sites });
  }
  rows.sort((a, b) => a.sites.length - b.sites.length || a.n.localeCompare(b.n));
  console.log(`${doc.name}: ${rows.length} of ${list.length} enzymes cut\n`);
  const unique = rows.filter((r) => r.sites.length === 1);
  if (unique.length) {
    console.log(`Unique cutters (${unique.length})`);
    for (const r of unique) {
      const e = ENZYME_BY_NAME[r.n];
      console.log(`  ${pad(r.n, 10)} ${pad(e.site, 14)} cuts after ${num(r.sites[0].top)}  ${overhangLabel(e)}`);
    }
  }
  const multi = rows.filter((r) => r.sites.length > 1);
  if (multi.length) {
    console.log(`\nMultiple cutters (${multi.length})`);
    for (const r of multi) console.log(`  ${pad(r.n, 10)} ${pad(r.sites.length + '×', 5)} at ${r.sites.map((s) => num(s.top)).join(', ')}`);
  }
  if (rest[1]) {
    const chosen = rest[1].split(',').map((s) => s.trim());
    const cuts = chosen.flatMap((n) => findSites(doc.seq, doc.circular, ENZYME_BY_NAME[n]));
    const frags = digest(doc.seq, doc.circular, cuts);
    console.log(`\nDigest with ${chosen.join(' + ')}: ${frags.length} fragment(s)`);
    for (const f of [...frags].sort((a, b) => b.len - a.len)) {
      const ends = `${f.left.end.type === 'blunt' ? 'blunt' : f.left.end.oh || 'end'} / ${f.right.end.type === 'blunt' ? 'blunt' : f.right.end.oh || 'end'}`;
      console.log(`  ${pad(num(f.len) + ' bp', 12)} ${num(f.start + 1)}..${num(f.end)}   ends ${ends}`);
    }
  }
}

function saveDoc(doc, out) {
  const ext = extname(out).toLowerCase();
  if (ext === '.dna') writeFileSync(out, writeDna(doc, { bindingSitesFor: (p) => findBindingSites(doc.seq, doc.circular, p.seq) }));
  else if (['.gb', '.gbk', '.genbank'].includes(ext)) writeFileSync(out, writeGenBank(doc));
  else if (['.fasta', '.fa', '.fna'].includes(ext)) writeFileSync(out, writeFasta(doc));
  else throw new Error(`Don't know how to write ${out}`);
}

function printClone(result, title) {
  const { doc, warnings } = result;
  const c = doc.cloning;
  console.log(title);
  console.log(
    `  fragments: ${c.inputs
      .map((i) => `${i.name} [${i.role}${i.orientation === 'reverse' ? ', reverse' : ''}] ${num(i.seqLen)} bp`)
      .join('  +  ')}`,
  );
  for (let n = 0; n < c.junctions.length; n++) {
    const j = c.junctions[n];
    const kind = j.type === 'overlap' ? `${j.overhang.length} bp overlap` : j.type === 'blunt' ? 'blunt' : `${j.type}′ ${j.overhang}`;
    const enz = j.enzymes.filter(Boolean).join('/') || '—';
    const sites = j.regeneratedSites.map((s) => `${s.name}@${num(s.pos + 1)}`).join(', ');
    const where = doc.circular && j.position === doc.seq.length ? 1 : j.position + 1;
    console.log(
      `  junction ${String(n + 1).padEnd(2)} ${pad(kind, 18)} ${pad(enz, 12)} at ${num(where)}` +
        `${sites ? `   regenerates ${sites}` : ''}`,
    );
  }
  console.log(`  product:   ${num(doc.seq.length)} bp ${doc.circular ? 'circular' : 'linear'}, ${doc.features.length} feature(s)`);
  for (const w of warnings) console.log(`  note: ${w}`);
}

try {
  if (!cmd || cmd === '-h' || cmd === '--help') {
    const src = readFileSync(new URL(import.meta.url)).toString().split('\n');
    const head = [];
    for (let i = 1; i < src.length && src[i].startsWith('//'); i++) head.push(src[i].replace(/^\/\/ ?/, ''));
    console.log(head.join('\n'));
    process.exit(0);
  }
  if (!file) throw new Error('Give a file to work on');
  const doc = load(file);
  const L = doc.seq.length;

  if (cmd === 'info') info(doc);
  else if (cmd === 'seq') {
    if (rest[0]) {
      const r = parseRegion(rest[0], L);
      const s = subseq(doc.seq, r.start, r.end, doc.circular);
      console.log(`${r.start + 1}..${r.end}  ${num(s.length)} bp  GC ${gcContent(s).toFixed(1)}%`);
      for (let i = 0; i < s.length; i += 60) console.log(`${String(r.start + i + 1).padStart(8)}  ${s.slice(i, i + 60).toLowerCase().match(/.{1,10}/g).join(' ')}`);
    } else {
      for (let i = 0; i < L; i += 60) console.log(`${String(i + 1).padStart(8)}  ${doc.seq.slice(i, i + 60).toLowerCase().match(/.{1,10}/g).join(' ')}`);
    }
  } else if (cmd === 'cut') cut(doc, rest[0]);
  else if (cmd === 'orfs') {
    const min = Number(rest[0] ?? 100);
    const orfs = findOrfs(doc.seq, doc.circular, min);
    console.log(`${orfs.length} ORF(s) of at least ${min} codons`);
    for (const o of orfs) console.log(`  ${pad(o.codons + ' codons', 14)} ${o.strand === 1 ? '+' : '-'} ${num(o.start + 1)}..${num(o.end)}`);
  } else if (cmd === 'primers') {
    const r = parseRegion(rest[0], L);
    const d = designPrimers(doc.seq, doc.circular, r.start, r.end, { targetTm: Number(rest[1] ?? 60) });
    for (const p of [d.forward, d.reverse]) {
      console.log(`${p.dir === 1 ? 'forward' : 'reverse'}  ${p.seq}`);
      console.log(`         ${p.seq.length} nt, Tm ${p.tm.toFixed(1)}°C, GC ${p.gc.toFixed(0)}%`);
    }
    console.log(`product  ${num(d.productLength)} bp`);
  } else if (cmd === 'pcr') {
    const { doc: product, warnings } = simulatePCR(doc, rest[0], rest[1]);
    for (const w of warnings) console.log(`note: ${w}`);
    info(product);
  } else if (cmd === 'clone') {
    const insertPath = rest[0];
    const enzymeSpec = rest[1];
    const out = rest[2];
    if (!insertPath || !enzymeSpec) {
      throw new Error('Usage: bg clone <vector> <insert> <enzyme,...> [out.dna|.gb|.fasta]');
    }
    const insertDoc = load(insertPath);
    const names = enzymeSpec.split(',').map((s) => s.trim()).filter(Boolean);
    const unknown = names.filter((n) => !ENZYME_BY_NAME[n]);
    if (unknown.length) throw new Error(`Unknown enzyme(s): ${unknown.join(', ')}`);
    const fragsOf = (d) =>
      digest(d.seq, d.circular, names.flatMap((n) => findSites(d.seq, d.circular, ENZYME_BY_NAME[n])));
    const vecFrags = fragsOf(doc);
    const insFrags = fragsOf(insertDoc);
    if (vecFrags.length < 2) {
      throw new Error(`${doc.name} is not cut into a backbone by ${names.join(', ')} (${vecFrags.length} fragment(s))`);
    }
    const backbone = [...vecFrags].sort((a, b) => b.len - a.len)[0];
    const cutEnds = (f) => f.left.end.type !== 'natural' && f.right.end.type !== 'natural';
    const insertFrag =
      !insertDoc.circular && insFrags.length >= 3
        ? insFrags[1]
        : [...insFrags].filter(cutEnds).sort((a, b) => b.len - a.len)[0];
    if (!insertFrag) throw new Error(`Could not pick an insert from ${insertDoc.name}; cut it so both ends are enzyme ends`);
    const result = restrictionClone({
      backbone: { doc, fragment: backbone },
      insert: { doc: insertDoc, fragment: insertFrag },
      circular: doc.circular,
    });
    printClone(result, `Restriction cloning with ${names.join(' + ')}`);
    if (out) {
      saveDoc(result.doc, out);
      console.log(`Wrote ${out}`);
    }
  } else if (cmd === 'gibson') {
    const args = rest.slice();
    const linear = args.includes('--linear');
    const withoutFlags = args.filter((a) => a !== '--linear');
    let out = null;
    const oi = withoutFlags.findIndex((a) => a === '-o' || a === '--out');
    if (oi >= 0) {
      out = withoutFlags[oi + 1];
      withoutFlags.splice(oi, 2);
    }
    if (!withoutFlags.length) {
      throw new Error('Usage: bg gibson <fragment1> <fragment2> [fragment3 ...] [-o out] [--linear]');
    }
    const fragments = [doc, ...withoutFlags.map((p) => load(p))];
    const result = gibsonAssemble(fragments, { circular: !linear, name: 'Gibson construct' });
    printClone(result, `Gibson assembly${linear ? ' (linear)' : ''}`);
    if (out) {
      saveDoc(result.doc, out);
      console.log(`Wrote ${out}`);
    }
  } else if (cmd === 'convert') {
    const out = rest[0];
    if (!out) throw new Error('Give an output file name');
    saveDoc(doc, out);
    console.log(`Wrote ${out}`);
  } else throw new Error(`Unknown command "${cmd}". Run with --help.`);
} catch (e) {
  console.error(`bg: ${e.message}`);
  process.exit(1);
}
