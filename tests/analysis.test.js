// Multiple alignment, protein properties, PDB, PDF, tracks, methylation, history
// jumps, collections' records and the export helpers.

import test from 'node:test';
import assert from 'node:assert/strict';
import { JSDOM } from 'jsdom';
import { newDocument } from '../js/model.js';
import { History } from '../js/undo.js';
import { multipleAlign, alignmentVariants } from '../js/msa.js';
import { proteinProperties, hydropathyProfile, hydrophobicSegments, proteinMotifs, chargeAt, isoelectricPoint } from '../js/protein.js';
import { parsePdb, chainSequence, normalise } from '../js/pdb.js';
import { imagePdf } from '../js/pdf.js';
import { gcProfile, gcSkewProfile, cumulativeSkew, summarize } from '../js/tracks.js';
import { methylationAt, methylationSummary, methylationKinds } from '../js/methylation.js';
import { ENZYME_BY_NAME, findSites } from '../js/enzymes.js';
import { recordFor } from '../js/collection.js';

// ---- multiple alignment ----------------------------------------------------------

test('every sequence survives alignment unchanged once the gaps are removed', () => {
  const seqs = ['ACGTACGTAAGGCCTTAAGGCCTTACGT', 'ACGTACGTAAGGCTTAAGGCCTTACGT', 'ACGTACCTAAGGCCTTAAGGGCCTTACGT', 'ACGTACGTAAGGCCTTAAGGCCTTAC'];
  const r = multipleAlign(seqs.map((seq, i) => ({ name: `s${i}`, seq })));
  assert.equal(new Set(r.rows.map((x) => x.length)).size, 1, 'all rows are the same width');
  r.rows.forEach((row, i) => assert.equal(row.replace(/-/g, ''), seqs[i]));
  assert.deepEqual(r.names, ['s0', 's1', 's2', 's3']);
  assert.equal(r.consensus.length, r.rows[0].length);
});

test('differences between aligned rows are classified and positioned along the first row', () => {
  const r = multipleAlign([{ name: 'ref', seq: 'ACGTACGTAAGGCCTTAAGGCCTTACGT' }, { name: 'del', seq: 'ACGTACGTAAGGCTTAAGGCCTTACGT' }, { name: 'sub', seq: 'ACGTACCTAAGGCCTTAAGGCCTTACGT' }]);
  const kinds = alignmentVariants(r.rows).map((v) => v.type).sort();
  assert.deepEqual(kinds, ['deletion', 'substitution']);
  const sub = alignmentVariants(r.rows).find((v) => v.type === 'substitution');
  assert.equal(sub.position, 6);
  assert.deepEqual(sub.bases, ['G', 'G', 'C']);
});

test('multiple alignment refuses what it cannot do', () => {
  assert.throws(() => multipleAlign([{ name: 'a', seq: 'ACGT' }]), /at least two/);
  assert.throws(() => multipleAlign([{ name: 'a', seq: 'A'.repeat(50) }, { name: 'b', seq: 'ACGT' }], { maxLength: 10 }), /limited/);
});

// ---- protein -----------------------------------------------------------------------

test('protein properties are in the right range and respond to composition', () => {
  const basic = proteinProperties('MKKRRKKRRKKR');
  const acidic = proteinProperties('MDDEEDDEEDDE');
  assert.ok(basic.pI > 10.5 && acidic.pI < 4.5);
  assert.ok(basic.chargeAt7 > 5 && acidic.chargeAt7 < -5);
  assert.ok(Math.abs(chargeAt('ACDE', isoelectricPoint('ACDE'))) < 0.01, 'charge is zero at the pI');
  const p = proteinProperties('GAVLIGAVLIGAVLI*TAIL');
  assert.equal(p.length, 15, 'translation stops at the first stop');
  assert.ok(p.gravy > 2, 'a Gly/Ala/Val/Leu/Ile protein is hydrophobic');
  assert.equal(proteinProperties('W').extinction, 5500);
  assert.equal(proteinProperties('').mass, 0);
});

test('hydropathy finds a membrane-like stretch, and motifs skip a site inside a tag', () => {
  const p = `MKKDD${'LLLIVVAALLIVVLLAAIIV'}DDEEKK`;
  const segs = hydrophobicSegments(p);
  assert.equal(segs.length, 1);
  assert.ok(segs[0].start >= 3 && segs[0].end <= 28);
  assert.ok(hydropathyProfile(p, 9).length > 0);
  const names = proteinMotifs('MDYKDDDDKGNSTAAHHHHHH').map((m) => m.name);
  assert.ok(names.includes('FLAG tag') && names.includes('His tag') && names.includes('N-glycosylation sequon'));
  assert.ok(!names.includes('Enterokinase site'), 'the site inside FLAG is part of it');
  assert.equal(proteinMotifs('NPS').some((m) => m.name.includes('glycosylation')), false, 'N-P-S is not a sequon');
});

// ---- PDB ---------------------------------------------------------------------------

const atom = (i, res, chain, n, x, y, z) => `ATOM  ${String(i).padStart(5)}  CA  ${res} ${chain}${String(n).padStart(4)}    ${x.toFixed(3).padStart(8)}${y.toFixed(3).padStart(8)}${z.toFixed(3).padStart(8)}  1.00  0.00           C`;

test('a PDB file gives its Cα trace, chains and declared secondary structure', () => {
  const text = ['TITLE     TEST', 'HELIX    1   1 ALA A    1  GLY A    2  1', 'SHEET    1   A 1 SER A   3  SER A   3  0',
    atom(1, 'ALA', 'A', 1, 0, 0, 0), atom(2, 'GLY', 'A', 2, 3.8, 0, 0), atom(3, 'SER', 'A', 3, 7.6, 0, 0), atom(4, 'LYS', 'B', 1, 0, 5, 0),
    'ATOM      5  N   LYS B   1       0.000   5.000   1.000  1.00  0.00           N', 'ENDMDL', atom(9, 'ALA', 'A', 9, 99, 99, 99)].join('\n');
  const s = parsePdb(text, { name: 'x' });
  assert.equal(s.title, 'TEST');
  assert.equal(s.chains.length, 2);
  assert.deepEqual(s.chains[0].atoms.map((a) => a.ss), ['helix', 'helix', 'sheet']);
  assert.equal(s.residues, 4, 'only Cα atoms, and nothing after ENDMDL');
  assert.equal(chainSequence(s.chains[0]), 'AGS');
  const { centre, radius } = normalise(s);
  assert.ok(radius > 3 && centre.length === 3);
  assert.throws(() => parsePdb('hello'), /No Cα/);
});

// ---- PDF ---------------------------------------------------------------------------

test('the PDF writer produces a one-page file whose cross-reference table is exact', () => {
  const jpeg = new Uint8Array([0xff, 0xd8, 0xff, 0xd9]);
  const bytes = imagePdf(jpeg, { width: 400, height: 300, title: 'Map (a) \\ b', subtitle: '510 bp' });
  const text = Buffer.from(bytes).toString('latin1');
  assert.ok(text.startsWith('%PDF-1.4'));
  assert.ok(text.trimEnd().endsWith('%%EOF'));
  const xref = Number(/startxref\n(\d+)/.exec(text)[1]);
  assert.equal(text.slice(xref, xref + 4), 'xref');
  const entries = [...text.slice(xref).matchAll(/(\d{10}) 00000 n /g)].map((m) => Number(m[1]));
  assert.equal(entries.length, 6);
  entries.forEach((offset, i) => assert.equal(text.slice(offset, offset + 7), `${i + 1} 0 obj`, `object ${i + 1} sits where the table says`));
  assert.match(text, /\/Count 1/);
  assert.match(text, /Map \\\(a\\\) \\\\ b/, 'parentheses and backslashes in the title are escaped');
  const wide = imagePdf(jpeg, { width: 900, height: 300 });
  assert.match(Buffer.from(wide).toString('latin1'), /MediaBox \[0 0 842 595\]/, 'a wide image gets a landscape page');
});

// ---- tracks --------------------------------------------------------------------------

test('GC content and skew follow the sequence, and circular windows wrap', () => {
  const seq = `${'GC'.repeat(50)}${'AT'.repeat(50)}`;
  const gc = gcProfile(seq, 20, 10);
  assert.equal(gc[1].value, 100);
  assert.equal(gc[gc.length - 2].value, 0);
  const edge = gcProfile(`GGGG${'AAAA'.repeat(10)}`, 8, 4, true);
  assert.equal(edge[0].value, 50, 'a window at the origin of a circle also sees the end');
  assert.equal(gcProfile(`GGGG${'AAAA'.repeat(10)}`, 8, 4, false)[0].value, 100, 'a linear edge window uses only the bases that exist');
  assert.equal(gcSkewProfile('GGGGGGGG', 4, 4)[0].value, 1);
  assert.equal(gcSkewProfile('CCCCCCCC', 4, 4)[0].value, -1);
  assert.equal(cumulativeSkew('GGCC', 1).map((p) => p.value).join(), '1,2,1,0');
  assert.deepEqual(summarize([{ value: 1 }, { value: 3 }]), { min: 1, max: 3, mean: 2 });
  assert.deepEqual(gcProfile('', 10, 1), []);
});

// ---- methylation -----------------------------------------------------------------------

test('methylation is judged from the sequence around each site', () => {
  const dam = (seq) => { const d = { seq, circular: false }; return methylationSummary(d, 'XbaI', findSites(seq, false, ENZYME_BY_NAME.XbaI)); };
  assert.deepEqual(dam('AAAATCTAGAAAAA'), [], 'XbaI alone is not affected');
  assert.deepEqual(dam('AAAATCTAGATCAAAA'), [{ kind: 'dam', effect: 'overlap' }], 'TCTAGATC overlaps GATC');
  assert.deepEqual(dam('AAGATCTAGAAA'), [{ kind: 'dam', effect: 'overlap' }], 'GATCTAGA overlaps on the other side');
  assert.deepEqual(dam('GATCTCTAGAAAA'), [], 'GATC next to the site, not overlapping it');
  assert.deepEqual(dam('GATCAAAATCTAGAAAAA'), [], 'a GATC further away does not matter');
  assert.deepEqual(methylationAt({ seq: 'AAGGATCCAA' }, 'BclI', { start: 2, end: 8 }).map((x) => x.effect), ['blocked']);
  assert.ok(methylationKinds('NotI').includes('cpg') && methylationKinds('BamHI').length === 0);
  const wrap = { seq: 'TCAAAAAAAATCTAGA', circular: true }; // the site ends in GA, and the sequence begins TC: GATC across the origin
  assert.equal(methylationSummary(wrap, 'XbaI', findSites(wrap.seq, true, ENZYME_BY_NAME.XbaI)).length, 1);
  assert.ok(!ENZYME_BY_NAME.BamHI.methylation.includes('dam') && ENZYME_BY_NAME.ClaI.methylation.includes('dam'));
});

// ---- undo timeline ------------------------------------------------------------------------

test('the undo history can be listed and jumped through in both directions', () => {
  const h = new History();
  const states = ['s0', 's1', 's2', 's3'];
  h.record(states[0], 'first');
  h.record(states[1], 'second');
  h.record(states[2], 'third');
  assert.deepEqual(h.timeline(), { labels: ['first', 'second', 'third'], position: 3 });
  assert.equal(h.jumpTo(1, 's3'), 's1', 'back to just after the first edit');
  assert.deepEqual(h.timeline(), { labels: ['first', 'second', 'third'], position: 1 });
  assert.equal(h.jumpTo(3, 's1'), 's3', 'and forward again');
  assert.equal(h.jumpTo(0, 's3'), 's0');
  assert.equal(h.jumpTo(99, 's0'), 's3', 'cannot go past the latest');
  h.jumpTo(1, 's3');
  h.record('s1', 'new', {});
  assert.deepEqual(h.timeline().labels, ['first', 'new'], 'a new edit drops the steps that were ahead');
});

// ---- collections ----------------------------------------------------------------------------

test('a saved record carries what search needs, and a copy of the file', () => {
  const doc = newDocument({ name: 'pXY', seq: 'ACGT'.repeat(10), circular: true });
  doc.features.push({ id: 'f', name: 'lacZ', type: 'CDS', strand: 1, color: '#000', segments: [{ start: 0, end: 12 }], qualifiers: [], extra: {} });
  doc.primers.push({ id: 'p', name: 'M13 fwd', seq: 'GTAAAACGACGGCCAGT', description: '', extra: {} });
  const r = recordFor(doc, '  Cloning  ');
  assert.equal(r.collection, 'Cloning');
  assert.deepEqual([r.featureNames, r.primerNames, r.bp, r.circular], [['lacZ'], ['M13 fwd'], 40, true]);
  r.doc.seq = 'X';
  assert.equal(doc.seq.length, 40, 'the record is a copy');
  assert.equal(recordFor(doc, '').collection, 'My files');
});

// ---- export helpers ------------------------------------------------------------------------------

test('exported SVG stands alone: no CSS variables, a white background, a size', async () => {
  const dom = new JSDOM('<!doctype html><body></body>');
  Object.assign(globalThis, { document: dom.window.document, XMLSerializer: dom.window.XMLSerializer, Node: dom.window.Node, HTMLElement: dom.window.HTMLElement });
  const { exportableSvg, printSheet } = await import('../js/ui/export.js');
  const svg = dom.window.document.createElementNS('http://www.w3.org/2000/svg', 'svg');
  svg.setAttribute('viewBox', '0 0 300 200');
  const rect = dom.window.document.createElementNS('http://www.w3.org/2000/svg', 'rect');
  rect.setAttribute('fill', 'var(--accent)');
  rect.setAttribute('stroke', 'var(--line)');
  rect.setAttribute('style', 'cursor:pointer;fill:var(--muted)');
  svg.append(rect);
  const { text, width, height } = exportableSvg(svg);
  assert.deepEqual([width, height], [300, 200]);
  assert.ok(!text.includes('var('), 'every variable is resolved');
  assert.match(text, /fill="#0d6e66"/);
  assert.match(text, /fill="#ffffff"/, 'white background');
  assert.match(text, /^<\?xml/);

  const doc = newDocument({ name: 'sheet', seq: 'ACGT'.repeat(40), circular: true });
  doc.features.push({ id: 'f', name: 'gene', type: 'CDS', strand: 1, color: '#e0884c', segments: [{ start: 10, end: 50 }], qualifiers: [], extra: {} });
  const sheet = printSheet(doc);
  assert.equal(sheet.querySelector('h1').textContent, 'sheet');
  assert.ok(sheet.querySelector('.printmap svg'));
  assert.equal(sheet.querySelector('pre.seq').textContent.split('\n').length, 3, '160 bp in lines of 60');
  assert.match(sheet.textContent, /gene/);
});
