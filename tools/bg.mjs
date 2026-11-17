#!/usr/bin/env node
// BublaGene command line: a way to use the library before there is a UI.
//
//   bg info      <file>                  summary, features, primers
//   bg seq       <file> [start..end]     sequence (1-based, inclusive)
//   bg cut       <file> [enzyme,...]     restriction sites and fragment sizes
//   bg orfs      <file> [minCodons]      open reading frames
//   bg primers   <file> <start..end>     design a primer pair for a region
//   bg pcr       <file> <fwd> <rev>      simulate PCR, print the product
//   bg design    <file> <start..end>     design with Tm/product/tail/restriction options
//   bg gel       <bp,bp,...>             fragment sizes on a simulated gel
//   bg enzymes   list|export|import ...  curated/custom enzyme sets
//   bg digest    <file> <enzyme,...>     per-enzyme fragment tables
//   bg overhangs <file> [enzymes]        Type IIS overhangs and compatibility
//   bg clone     <vector> <insert> <enzyme,...> [out]   restriction cloning
//   bg gibson    <f1> <f2> [...] [-o out] [--linear]    Gibson assembly
//   bg assemble  <f1> <f2> ...          N-fragment homology assembly (--plan FILE)
//   bg golden-gate <f1> <f2> ...        Type IIS assembly (--enzyme BsaI, --plan FILE)
//   bg assembly-primers <f1> <f2> ...   overlap primers for each fragment
//   bg plan      <file>                 print and verify an assembly/cloning plan
//   bg replay    <plan.json>            rebuild a product from its plan
//   bg batch     <library.json> <designs.json>   many constructs in one run
//   bg align     <a> <b> [--local] [--band N]   pairwise alignment
//   bg dotplot   <a> <b>                        dot-plot grid
//   bg contig    <r1> <r2> ...                  assemble reads into contigs
//   bg variants  <ref> <query> [--features]     SNPs/indels between two sequences
//   bg trace     <trace.ab1|.scf>               chromatogram summary
//   bg trim      <trace> [--mode mott|phred] [--cutoff 20] [--out out.fasta]
//   bg history   <file>                  read the cloning record beside a file
//   bg translate <file> [start..end] [table]  translate with an NCBI code
//   bg codons                              list and validate genetic codes
//   bg annotate  <file> [out]              recompute translations and derived qualifiers
//   bg composition <file>                  base composition and GC
//   bg convert   <file> <out.gb|.fasta|.dna>
//
// Files may be .dna (binary), .gb/.gbk/.genbank, or .fasta/.fa.

import { readFileSync, writeFileSync, existsSync } from 'node:fs';
import { basename, extname } from 'node:path';
import { readDna, writeDna, describeSegments } from '../js/formats/dna.js';
import { readGenBank, writeGenBank, readFasta, writeFasta } from '../js/formats/genbank.js';
import { ENZYMES, ENZYME_BY_NAME, COMMON_SET, findSites, digest, cutEnd, overhangLabel, enzymeInfo, enzymesToJSON, enzymesFromJSON, mergeEnzymeSets, uniqueCutters, cuttersByCount, fragmentTable, typeIISOverhangs, overhangMatrix } from '../js/enzymes.js';
import { findBindingSites, designPrimers, designPrimerPairs, simulatePCR, pcrDiagnostics, gelView } from '../js/primers.js';
import { restrictionClone, gibsonAssemble, cloningRecord, cloningFromJSON, cloningSidecarName } from '../js/cloning.js';
import { homologyAssemble, goldenGateAssemble, designAssemblyPrimers, assemblyPlan, assemblyFromJSON, replayAssembly, batchAssembly } from '../js/assembly.js';
import { alignGlobal, alignLocal, alignPair, alignBanded, dotPlot, dotPlotText, formatAlignment } from '../js/align.js';
import { assembleContigs } from '../js/contig.js';
import { callVariants, variantsToFeatures, formatVariant } from '../js/variants.js';
import { readTrace, trimTrace, qualityStats, traceToDocument } from '../js/traces.js';
import { codingSequence, featureLength, featureSpan, isCodingSegment, translateFeature, featureTableId, featureProteinMass, recomputeFeatures, qualifierText } from '../js/model.js';
import { gcContent, meltingTemp, translate, findOrfs, subseq, rangeLength, nucleotideCount } from '../js/seq.js';
import { TABLE_IDS, CODON_TABLES, validateAllTables } from '../js/codon-tables.js';

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

// Split "--flag value" pairs from positional arguments.
function splitFlags(args, known) {
  const positional = [];
  const flags = {};
  for (let i = 0; i < args.length; i++) {
    const a = args[i];
    if (a.startsWith('--')) {
      const key = a.slice(2);
      if (known.includes(key)) {
        const next = args[i + 1];
        if (next !== undefined && !next.startsWith('--')) { flags[key] = args[++i]; } else { flags[key] = true; }
        continue;
      }
      flags[key] = true;
      continue;
    }
    positional.push(a);
  }
  return { positional, flags };
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
        const prot = translateFeature(doc, f);
        const table = featureTableId(f, doc);
        const flags = [];
        if (cds.length % 3) flags.push('not a whole number of codons');
        if (prot.slice(0, -1).includes('*')) flags.push('internal stop');
        if (table !== 1) flags.push(`code ${table}`);
        const aa = prot.length - (prot.endsWith('*') ? 1 : 0);
        extra = `  ${aa} aa, ${Math.round(featureProteinMass(doc, f)).toLocaleString('en-US')} Da${flags.length ? ` (${flags.join(', ')})` : ''}`;
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

function cut(doc, names, flags = {}) {
  let list = (names ? names.split(',') : COMMON_SET).map((n) => n.trim()).filter(Boolean);
  if (flags.n != null) list = cuttersByCount(doc.seq, doc.circular, list, Number(flags.n)).map((r) => r.enzyme.name);
  if (flags.unique) list = uniqueCutters(doc.seq, doc.circular, list).map((r) => r.enzyme.name);
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

function printRecord(record, title, size, warnings = [], features = null) {
  console.log(title);
  console.log(
    `  fragments: ${record.inputs
      .map((i) => `${i.name} [${i.role}${i.orientation === 'reverse' ? ', reverse' : ''}] ${num(i.seqLen)} bp`)
      .join('  +  ')}`,
  );
  for (let n = 0; n < record.junctions.length; n++) {
    const j = record.junctions[n];
    const kind = j.type === 'overlap' ? `${j.overhang.length} bp overlap` : j.type === 'blunt' ? 'blunt' : `${j.type}′ ${j.overhang}`;
    const enz = j.enzymes.filter(Boolean).join('/') || '—';
    const sites = j.regeneratedSites.map((s) => `${s.name}@${num(s.pos + 1)}`).join(', ');
    const where = record.circular && j.position === size ? 1 : j.position + 1;
    console.log(
      `  junction ${String(n + 1).padEnd(2)} ${pad(kind, 18)} ${pad(enz, 12)} at ${num(where)}` +
        `${sites ? `   regenerates ${sites}` : ''}`,
    );
  }
  console.log(`  product:   ${num(size)} bp ${record.circular ? 'circular' : 'linear'}${features == null ? '' : `, ${features} feature(s)`}`);
  for (const w of warnings) console.log(`  note: ${w}`);
}

function printClone(result, title) {
  printRecord(result.doc.cloning, title, result.doc.seq.length, result.warnings, result.doc.features.length);
}

function saveRecord(doc, out) {
  const record = cloningRecord(doc);
  if (!record) return null;
  const path = cloningSidecarName(out);
  writeFileSync(path, `${JSON.stringify(record, null, 2)}\n`);
  console.log(`Wrote ${path}`);
  return path;
}

try {
  if (!cmd || cmd === '-h' || cmd === '--help') {
    const src = readFileSync(new URL(import.meta.url)).toString().split('\n');
    const head = [];
    for (let i = 1; i < src.length && src[i].startsWith('//'); i++) head.push(src[i].replace(/^\/\/ ?/, ''));
    console.log(head.join('\n'));
    process.exit(0);
  }
  if (cmd === 'codons') {
    const checks = validateAllTables();
    console.log(`NCBI genetic codes available: ${TABLE_IDS.join(', ')}`);
    for (const v of checks) console.log(`  ${String(v.id).padStart(2)}  ${pad(v.name, 48)} ${v.ok ? 'ok' : `INVALID: ${v.issues.join(', ')}`}`);
    process.exit(0);
  }
  if (cmd === 'gel') {
    const { positional, flags } = splitFlags(rest, ['ladder', 'agarose']);
    const sizes = (file ?? positional[0] ?? '').split(',').map((n) => Number(n.trim())).filter((n) => n > 0);
    if (!sizes.length) throw new Error('Give fragment sizes in bp, comma-separated, e.g. bg gel 1000,500,200');
    const gel = gelView(sizes, { ladder: flags.ladder ?? '1 kb ladder', agarose: Number(flags.agarose ?? 1) });
    for (const b of gel.bands) console.log(`  ${pad(num(b.bp) + ' bp', 14)} ${(b.migration * 100).toFixed(0).padStart(3)}% down the gel`);
    console.log(`\n${gel.bands.length} band(s); ladder ${gel.ladder}`);
    process.exit(0);
  }
  if (cmd === 'enzymes') {
    const sub = file;
    const { positional, flags } = splitFlags(rest, ['file', 'type', 'search', 'json', 'merge', 'out']);
    const custom = flags.file ? enzymesFromJSON(readFileSync(flags.file, 'utf8')) : [];
    if (sub === 'list' || !sub) {
      let list = custom.length ? mergeEnzymeSets(ENZYMES, custom) : ENZYMES;
      if (flags.type) list = list.filter((e) => e.type === String(flags.type).toUpperCase());
      if (flags.search) list = list.filter((e) => e.name.toLowerCase().includes(String(flags.search).toLowerCase()) || e.site.includes(String(flags.search).toUpperCase()));
      if (flags.json) {
        console.log(JSON.stringify(enzymesToJSON(list), null, 2));
      } else {
        for (const e of list) {
          console.log(`  ${pad(e.name, 11)} ${pad(e.site, 14)} ${pad(`cut ${e.top}/${e.bottom}`, 12)} ${pad(e.type, 4)} ${pad(e.overhang === 'blunt' ? 'blunt' : `${e.overhang}′ ${e.overhangLength}`, 8)} ${e.methylation.length ? `methylation: ${e.methylation.join(',')}` : ''}`);
        }
        console.log(`\n${list.length} enzyme(s)${custom.length ? ` (${custom.length} from ${flags.file})` : ''}`);
      }
    } else if (sub === 'export') {
      const out = positional[0];
      if (!out) throw new Error('Usage: bg enzymes export <out.json> [name,...]');
      const names = positional[1] ? positional[1].split(',').map((x) => x.trim()).filter(Boolean) : null;
      const list = names ? names.map((n) => {
        const e = ENZYME_BY_NAME[n] ?? custom.find((x) => x.name === n);
        if (!e) throw new Error(`Unknown enzyme ${n}`);
        return e;
      }) : (custom.length ? mergeEnzymeSets(ENZYMES, custom) : ENZYMES);
      writeFileSync(out, `${JSON.stringify(enzymesToJSON(list), null, 2)}\n`);
      console.log(`Wrote ${out} (${list.length} enzyme(s))`);
    } else if (sub === 'import') {
      const setPath = positional[0];
      if (!setPath) throw new Error('Usage: bg enzymes import <set.json> [--out merged.json] [--merge base.json]');
      const incoming = enzymesFromJSON(readFileSync(setPath, 'utf8'));
      let base = ENZYMES;
      if (flags.merge) base = enzymesFromJSON(readFileSync(flags.merge, 'utf8'));
      const merged = mergeEnzymeSets(base, incoming);
      if (flags.out) {
        writeFileSync(flags.out, `${JSON.stringify(enzymesToJSON(merged), null, 2)}\n`);
        console.log(`Wrote ${flags.out} (${merged.length} enzyme(s))`);
      } else {
        for (const e of incoming) console.log(`  ${pad(e.name, 11)} ${pad(e.site, 12)} cut ${e.top}/${e.bottom}  ${e.type}`);
        console.log(`\n${incoming.length} enzyme(s) validated; ${merged.length} after merging.`);
      }
    } else {
      throw new Error(`Unknown enzymes subcommand "${sub}". Use list, export or import.`);
    }
    process.exit(0);
  }
  if (cmd === 'replay') {
    if (!file) throw new Error('Usage: bg replay <plan.json>');
    const plan = assemblyFromJSON(readFileSync(file, 'utf8'));
    const out = replayAssembly(plan);
    console.log(`${out.ok ? 'OK' : 'MISMATCH'}: rebuilt ${out.seq.length} bp, digest ${out.digest} (expected ${plan.product.digest})`);
    if (!out.ok) process.exit(1);
    process.exit(0);
  }
  if (cmd === 'batch') {
    if (!file) throw new Error('Usage: bg batch <library.json> <designs.json> [--out-dir DIR]');
    const { positional, flags } = splitFlags(rest, ['out-dir']);
    const designsPath = positional[0];
    if (!designsPath) throw new Error('Give a designs JSON file');
    const rawLib = JSON.parse(readFileSync(file, 'utf8'));
    const library = {};
    for (const [name, value] of Object.entries(rawLib)) {
      if (typeof value === 'string') library[name] = { name, seq: value };
      else if (value.file) library[name] = load(value.file);
      else library[name] = { seq: value.seq, name };
    }
    const designs = JSON.parse(readFileSync(designsPath, 'utf8'));
    const batch = batchAssembly(library, Array.isArray(designs) ? designs : designs.designs);
    if (flags['out-dir']) {
      const { mkdirSync } = await import('node:fs');
      mkdirSync(flags['out-dir'], { recursive: true });
      for (const r of batch.results.filter((x) => x.ok)) {
        const out = `${flags['out-dir']}/${r.name.replace(/[^\w.-]+/g, '_')}.gb`;
        writeFileSync(out, writeGenBank(r.doc));
        writeFileSync(`${flags['out-dir']}/${r.name.replace(/[^\w.-]+/g, '_')}.assembly.json`, `${JSON.stringify(r.plan, null, 2)}\n`);
      }
    }
    console.log(JSON.stringify({
      ok: batch.ok, failed: batch.failed,
      results: batch.results.map((r) => ({ name: r.name, ok: r.ok, length: r.doc?.seq.length, error: r.error })),
    }, null, 2));
    process.exit(batch.failed ? 1 : 0);
  }
  if (cmd === 'trace' || cmd === 'trim') {
    if (!file) throw new Error(`Usage: bg ${cmd} <trace.ab1|.scf> [flags]`);
    const { flags } = splitFlags(rest, ['mode', 'cutoff', 'window', 'out', 'json']);
    const trace = await readTrace(readFileSync(file), basename(file));
    const stats = qualityStats(trace);
    const trim = trimTrace(trace, { mode: flags.mode ?? 'mott', cutoff: Number(flags.cutoff ?? 20), window: Number(flags.window ?? 10) });
    if (flags.json) {
      console.log(JSON.stringify({ format: trace.format, name: trace.name, bases: trace.bases.length, samples: trace.sampleCount, quality: stats, trim }, null, 2));
    } else {
      console.log(`${trace.format} trace ${trace.name || basename(file)}: ${num(trace.bases.length)} bases, ${num(trace.sampleCount)} samples`);
      console.log(`  quality: mean ${stats.mean.toFixed(1)}, min ${stats.min}, max ${stats.max}`);
      console.log(`  trim (${trim.mode}, cutoff ${trim.cutoff}): keep ${trim.start + 1}..${trim.end} (${trim.length} bases), mean ${trim.mean.toFixed(1)}`);
      console.log(`  ${trace.bases.slice(trim.start, Math.min(trim.end, trim.start + 100))}${trim.length > 100 ? '…' : ''}`);
    }
    if (cmd === 'trim' && flags.out) {
      const doc = traceToDocument(trace, { trimmed: false });
      doc.seq = trace.bases.slice(trim.start, trim.end);
      writeFileSync(flags.out, writeFasta(doc));
      console.log(`Wrote ${flags.out} (${doc.seq.length} bases)`);
    }
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
  } else if (cmd === 'cut') {
    const { positional, flags } = splitFlags(rest, ['n', 'unique']);
    cut(doc, positional[0], flags);
  }
  else if (cmd === 'digest') {
    const { positional } = splitFlags(rest, []);
    const names = (positional[0] ?? '').split(',').map((x) => x.trim()).filter(Boolean);
    if (!names.length) throw new Error('Usage: bg digest <file> <enzyme,...>');
    const rows = fragmentTable(doc.seq, doc.circular, names);
    for (const row of rows) {
      console.log(`${row.name}  ${row.count} site(s)`);
      for (const f of row.fragments) {
        const ends = `${f.left.type === 'blunt' ? 'blunt' : f.left.oh || 'end'} / ${f.right.type === 'blunt' ? 'blunt' : f.right.oh || 'end'}`;
        console.log(`  ${pad(num(f.len) + ' bp', 12)} ${num(f.start + 1)}..${num(f.end)}   ends ${ends}`);
      }
      if (!row.fragments.length) console.log('  (no cuts)');
    }
  }
  else if (cmd === 'overhangs') {
    const { positional, flags } = splitFlags(rest, ['json']);
    const names = (positional[0] ? positional[0].split(',') : ['BsaI', 'BsmBI', 'BbsI', 'SapI']).map((x) => x.trim()).filter(Boolean);
    const matrix = overhangMatrix(doc.seq, doc.circular, names);
    if (flags.json) {
      console.log(JSON.stringify({ overhangs: matrix.overhangs, rows: matrix.rows, matrix: matrix.matrix }, null, 2));
    } else {
      for (const row of matrix.rows) {
        console.log(`  ${pad(row.enzyme, 9)} ${row.strand === 1 ? '+' : '-'} ${pad(num(row.pos + 1), 8)} ${row.type}′ ${row.overhang} (rc ${row.reverseComplement})`);
      }
      console.log(`\nOverhangs: ${matrix.overhangs.join(', ') || 'none'}`);
      for (const oh of matrix.overhangs) console.log(`  ${oh}: ${matrix.groups[oh].length} site(s)`);
    }
  }
  else if (cmd === 'translate') {
    const table = Number(rest[1] ?? 1);
    let s = doc.seq;
    let where = 'whole sequence';
    if (rest[0] && rest[0].includes('..')) {
      const r = parseRegion(rest[0], L);
      s = subseq(doc.seq, r.start, r.end, doc.circular);
      where = `${r.start + 1}..${r.end}`;
    }
    const p = translate(s, table);
    console.log(`${where}  ${p.length} aa  genetic code ${table}`);
    for (let i = 0; i < p.length; i += 60) console.log(`${String(i + 1).padStart(7)}  ${p.slice(i, i + 60)}`);
  } else if (cmd === 'annotate') {
    recomputeFeatures(doc);
    const out = rest[0];
    if (!out) {
      for (const f of doc.features) {
        const t = qualifierText(f, 'translation');
        if (t) console.log(`${pad(f.type, 14)} ${pad(f.name, 22)} ${t.length} aa  ${t.slice(0, 40)}${t.length > 40 ? '…' : ''}`);
      }
      console.log(`${doc.features.length} feature(s) recomputed (no output file given)`);
    } else {
      saveDoc(doc, out);
      console.log(`Wrote ${out}`);
    }
  } else if (cmd === 'composition') {
    const c = nucleotideCount(doc.seq);
    console.log(`${doc.name}: ${c.length} bp, GC ${c.gcPercent.toFixed(1)}%`);
    for (const k of ['A', 'C', 'G', 'T', 'U', 'N']) if (c[k]) console.log(`  ${k}: ${c[k]}`);
  }
  else if (cmd === 'orfs') {
    const min = Number(rest[0] ?? 100);
    const orfs = findOrfs(doc.seq, doc.circular, min);
    console.log(`${orfs.length} ORF(s) of at least ${min} codons`);
    for (const o of orfs) console.log(`  ${pad(o.codons + ' codons', 14)} ${o.strand === 1 ? '+' : '-'} ${num(o.start + 1)}..${num(o.end)}`);
  } else if (cmd === 'primers' || cmd === 'design') {
    const { positional, flags } = splitFlags(rest, ['tm', 'min-product', 'max-product', 'fwd-tail', 'rev-tail', 'fwd-enzyme', 'rev-enzyme', 'flex', 'min-len', 'max-len', 'json']);
    const r = parseRegion(positional[0], L);
    const targetTm = Number(flags.tm ?? positional[1] ?? 60);
    const options = {
      targetTm,
      minTm: targetTm - 5,
      maxTm: targetTm + 8,
      minProduct: flags['min-product'] ? Number(flags['min-product']) : 0,
      maxProduct: flags['max-product'] ? Number(flags['max-product']) : Infinity,
      fwdTail: flags['fwd-tail'] ?? '',
      revTail: flags['rev-tail'] ?? '',
      fwdRestriction: flags['fwd-enzyme'] ?? null,
      revRestriction: flags['rev-enzyme'] ?? null,
      flex: Number(flags.flex ?? 0),
      minLen: Number(flags['min-len'] ?? 18),
      maxLen: Number(flags['max-len'] ?? 32),
    };
    const d = designPrimerPairs(doc.seq, doc.circular, r.start, r.end, options);
    if (flags.json) {
      console.log(JSON.stringify({
        forward: { sequence: d.forward.seq, annealing: d.forward.annealing, tm: d.forward.tm, gc: d.forward.gc },
        reverse: { sequence: d.reverse.seq, annealing: d.reverse.annealing, tm: d.reverse.tm, gc: d.reverse.gc },
        productLength: d.productLength, tmDifference: d.tmDiff, warnings: d.warnings,
      }, null, 2));
    } else {
      for (const p of [d.forward, d.reverse]) {
        console.log(`${p.dir === 1 ? 'forward' : 'reverse'}  ${p.seq}`);
        const props = p.properties;
        console.log(`         ${p.seq.length} nt, Tm ${p.tm.toFixed(1)}°C, GC ${p.gc.toFixed(0)}%, 3′ ${props.clamp.threePrimeGc ? 'G/C clamp' : 'no clamp'}`);
        if (props.hairpin.dg < -3) console.log(`         hairpin ΔG37 ${props.hairpin.dg.toFixed(1)} kcal/mol (${props.hairpin.stem} bp stem)`);
        if (props.selfDimer.dg < -5) console.log(`         self-dimer ΔG37 ${props.selfDimer.dg.toFixed(1)} kcal/mol`);
      }
      console.log(`product  ${num(d.productLength)} bp  (Tm difference ${d.tmDiff.toFixed(1)} °C)`);
      for (const w of d.warnings) console.log(`note: ${w}`);
    }
  } else if (cmd === 'pcr') {
    const { positional, flags } = splitFlags(rest, ['mismatch', 'anneal', 'three-prime', 'json']);
    const [fwdSeq, revSeq] = positional;
    if (!fwdSeq || !revSeq) throw new Error('Usage: bg pcr <file> <forward> <reverse> [--mismatch N] [--anneal N] [--json]');
    const opts = {
      minAnneal: Number(flags.anneal ?? 15),
      maxMismatches: Number(flags.mismatch ?? 0),
      threePrimeExact: Number(flags['three-prime'] ?? 1),
    };
    try {
      const result = simulatePCR(doc, fwdSeq, revSeq, opts);
      if (flags.json) {
        console.log(JSON.stringify({
          product: { length: result.doc.seq.length, sequence: result.doc.seq },
          junctions: result.junctions, warnings: result.warnings,
        }, null, 2));
      } else {
        for (const w of result.warnings) console.log(`note: ${w}`);
        info(result.doc);
        const gel = gelView([{ len: result.doc.seq.length }]);
        console.log(`\ngel: ${gel.list.join(', ')} bp`);
      }
    } catch (e) {
      if (e.diagnostics) for (const r of e.diagnostics.reasons) console.log(`reason: ${r}`);
      throw e;
    }
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
      saveRecord(result.doc, out);
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
      saveRecord(result.doc, out);
    }
  } else if (cmd === 'align') {
    const { positional, flags } = splitFlags(rest, ['band', 'local', 'global', 'json', 'width']);
    const other = positional.find((p) => !p.startsWith('-'));
    if (!other) throw new Error('Usage: bg align <a> <b> [--local] [--band N]');
    const b = load(other).seq;
    const opts = {
      local: !!flags.local,
      band: flags.band != null ? Number(flags.band) : null,
    };
    const result = (opts.local ? alignLocal : alignPair)(doc.seq, b, opts);
    if (flags.json) console.log(JSON.stringify({ score: result.score, identity: result.identity, matches: result.matches, mismatches: result.mismatches, gaps: result.gaps, alignedA: result.alignedA, alignedB: result.alignedB }, null, 2));
    else console.log(formatAlignment(result, { width: Number(flags.width ?? 60) }));
  } else if (cmd === 'dotplot') {
    const { positional } = splitFlags(rest, []);
    const other = positional.find((p) => !p.startsWith('-'));
    if (!other) throw new Error('Usage: bg dotplot <a> <b>');
    const b = load(other).seq;
    console.log(dotPlotText(doc.seq, b));
    const plot = dotPlot(doc.seq, b);
    console.log(`${plot.points.length} point(s), ${plot.width}×${plot.height}`);
  } else if (cmd === 'contig') {
    const { positional, flags } = splitFlags(rest, ['min-overlap', 'json']);
    const reads = [doc, ...positional.filter((p) => !p.startsWith('-')).map((p) => load(p))];
    const contigs = assembleContigs(reads, { minOverlap: Number(flags['min-overlap'] ?? 30) });
    if (flags.json) console.log(JSON.stringify(contigs.map((c) => ({ length: c.length, depth: c.depth, reads: c.reads, seq: c.seq })), null, 2));
    else {
      for (const [n, c] of contigs.entries()) {
        console.log(`contig ${n + 1}: ${num(c.length)} bp, ${c.reads.length} read(s), mean depth ${c.depth}`);
        console.log(`  ${c.seq.slice(0, 120)}${c.seq.length > 120 ? '…' : ''}`);
      }
      console.log(`${contigs.length} contig(s)`);
    }
  } else if (cmd === 'variants') {
    const { positional, flags } = splitFlags(rest, ['circular', 'json', 'features', 'out']);
    const other = positional.find((p) => !p.startsWith('-'));
    if (!other) throw new Error('Usage: bg variants <reference> <query> [--circular] [--features out.gb]');
    const query = load(other);
    const result = callVariants(doc.seq, query.seq, { circular: !!flags.circular || doc.circular });
    if (flags.json) {
      console.log(JSON.stringify({ identity: result.identity, variants: result.variants.map((v) => ({ type: v.type, position: v.position, ref: v.ref, alt: v.alt })) }, null, 2));
    } else {
      console.log(`${result.variants.length} variant(s); identity ${result.identity.toFixed(1)}%`);
      for (const v of result.variants) console.log(`  ${formatVariant(v)}`);
    }
    if (flags.features) {
      const out = load(other);
      out.features = [...out.features, ...variantsToFeatures(out, result.variants)];
      saveDoc(out, flags.features);
      console.log(`Wrote ${flags.features}`);
    }
  }
  else if (cmd === 'assemble') {
    const { positional, flags } = splitFlags(rest, ['min-overlap', 'name', 'plan', 'json', 'linear', 'circular']);
    const out = positional.find((p) => p.endsWith('.dna') || p.endsWith('.gb') || p.endsWith('.fasta'));
    let paths = positional.filter((p) => p !== out);
    if (out) paths = paths.slice(0, -1);
    const fragments = [doc, ...paths.filter((p) => !p.startsWith('-')).map((p) => load(p))];
    const result = homologyAssemble(fragments, {
      circular: !flags.linear, minOverlap: Number(flags['min-overlap'] ?? 15),
      name: flags.name ?? 'Assembly',
    });
    printRecord(result.doc.cloning, `Homology assembly (${fragments.length} fragments)`, result.doc.seq.length, result.warnings, result.doc.features.length);
    if (flags.json) console.log(JSON.stringify(result.plan, null, 2));
    if (out) {
      saveDoc(result.doc, out);
      console.log(`Wrote ${out}`);
      saveRecord(result.doc, out);
    }
    if (flags.plan) { writeFileSync(flags.plan, `${JSON.stringify(result.plan, null, 2)}\n`); console.log(`Wrote ${flags.plan}`); }
  } else if (cmd === 'golden-gate' || cmd === 'goldengate') {
    const { positional, flags } = splitFlags(rest, ['enzyme', 'name', 'plan', 'json', 'linear']);
    const fragments = [doc, ...positional.filter((p) => !p.startsWith('-')).map((p) => load(p))];
    const result = goldenGateAssemble(fragments, {
      circular: !flags.linear, enzyme: flags.enzyme ?? 'BsaI', name: flags.name ?? 'Golden Gate construct',
    });
    for (const c of result.checks) console.log(`  ${c.severity}: ${c.message}`);
    printRecord(result.doc.cloning, `Golden Gate assembly with ${flags.enzyme ?? 'BsaI'}`, result.doc.seq.length, result.warnings, result.doc.features.length);
    if (flags.json) console.log(JSON.stringify(result.plan, null, 2));
    if (flags.plan) { writeFileSync(flags.plan, `${JSON.stringify(result.plan, null, 2)}\n`); console.log(`Wrote ${flags.plan}`); }
  } else if (cmd === 'assembly-primers') {
    const { positional, flags } = splitFlags(rest, ['overlap', 'tm']);
    const fragments = [doc, ...positional.filter((p) => !p.startsWith('-')).map((p) => load(p))];
    const primers = designAssemblyPrimers(fragments, { overlap: Number(flags.overlap ?? 25), targetTm: Number(flags.tm ?? 60) });
    for (const p of primers) {
      console.log(`${p.name}`);
      console.log(`  forward  ${p.forward.sequence}  (tail ${p.forward.tail || '—'}, Tm ${p.forward.tm.toFixed(1)}°C)`);
      console.log(`  reverse  ${p.reverse.sequence}  (tail ${p.reverse.tail || '—'}, Tm ${p.reverse.tm.toFixed(1)}°C)`);
    }
  } else if (cmd === 'plan') {
    const a = assemblyPlan(doc);
    const c = cloningRecord(doc);
    if (a) {
      console.log(`Assembly plan (${a.method}, ${a.fragments.length} fragments, ${a.product.seqLen} bp)`);
      for (const j of a.junctions) console.log(`  junction at ${j.position + 1}: ${j.overhang?.length ?? 0} bp ${j.type}, scar ${j.scar ?? '—'}`);
      console.log(`  product digest ${a.product.digest}`);
      console.log(`  replay: ${replayAssembly(a).ok ? 'identical' : 'MISMATCH'}`);
    } else if (c) {
      printRecord(c, `Cloning record for ${doc.name}`, doc.seq.length);
    } else {
      console.log('No assembly or cloning plan on this document.');
    }
  }
  else if (cmd === 'history') {
    const path = cloningSidecarName(file);
    if (!existsSync(path)) throw new Error(`No cloning record beside ${file} (looked for ${path})`);
    const record = cloningFromJSON(readFileSync(path, 'utf8'));
    printRecord(record, `Cloning record for ${basename(file)}`, record.product.seqLen);
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
