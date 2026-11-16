import test from 'node:test';
import assert from 'node:assert/strict';
import { revcomp, translate, gcContent, meltingTemp, findPattern, findOrfs, subseq, rangeLength } from '../js/seq.js';
import { newDocument, replaceRange, rotateOrigin, reverseComplementDoc, featureSequence, featuresInRegion, newId } from '../js/model.js';
import { ENZYME_BY_NAME, findSites, digest, cutEnd, migration } from '../js/enzymes.js';
import { findBindingSites, designPrimers, simulatePCR } from '../js/primers.js';
import { parseXML, toXML, el } from '../js/xml.js';
import { parseLocation, formatLocation, readGenBank, writeGenBank, readFasta, writeFasta } from '../js/formats/genbank.js';

const feat = (name, start, end, strand = 1) => ({ id: newId(), name, type: 'misc_feature', strand, color: '#888', segments: [{ start, end }], qualifiers: [], extra: {} });

test('reverse complement and IUPAC', () => {
  assert.equal(revcomp('ATGCN'), 'NGCAT');
  assert.equal(revcomp('RYKM'), 'KMRY');
  assert.equal(revcomp('acgT'), 'Acgt');
});

test('translation', () => {
  assert.equal(translate('ATGGCCTAA'), 'MA*');
  assert.equal(translate('ATGGC'), 'M');
});

test('GC and Tm are sensible', () => {
  assert.equal(gcContent('GGCC'), 100);
  assert.equal(gcContent('ATGC'), 50);
  const tm = meltingTemp('AGCGGATAACAATTTCACACAGGA');
  assert.ok(tm > 55 && tm < 70, `Tm ${tm}`);
  assert.ok(meltingTemp('GCGCGCGCGCGCGCGCGC') > meltingTemp('ATATATATATATATATAT'));
});

test('pattern search on both strands and across the origin', () => {
  const seq = 'TTCAAAAAAAGAATTCAAAAGGTCTCAAAAAGAA'; // EcoRI site wraps around the end
  const hits = findPattern(seq, 'GAATTC', true);
  assert.deepEqual(hits.map((h) => h.start).sort((a, b) => a - b), [10, 31]);
  const bsa = findPattern('AAAGAGACCAAA', 'GGTCTC', false);
  assert.equal(bsa.length, 1);
  assert.equal(bsa[0].strand, -1);
});

test('restriction sites, cut positions and overhangs', () => {
  const seq = 'AAAAGAATTCAAAACTGCAGAAAA';
  const eco = findSites(seq, false, ENZYME_BY_NAME.EcoRI);
  assert.equal(eco.length, 1);
  assert.equal(eco[0].top, 5);
  assert.deepEqual(cutEnd(seq, false, eco[0]), { type: '5', oh: 'AATT' });
  const pst = findSites(seq, false, ENZYME_BY_NAME.PstI);
  assert.deepEqual(cutEnd(seq, false, pst[0]), { type: '3', oh: 'TGCA' });
  const sma = findSites('AACCCGGGAA', false, ENZYME_BY_NAME.SmaI);
  assert.deepEqual(cutEnd('AACCCGGGAA', false, sma[0]), { type: 'blunt', oh: '' });
});

test('type IIS enzyme cuts outside its site, on either strand', () => {
  //            site GGTCTC at 2, top cut = 2 + 7 = 9, bottom = 13
  const seq = 'AAGGTCTCAGCTTAAAAAAAA';
  const fwd = findSites(seq, false, ENZYME_BY_NAME.BsaI);
  assert.equal(fwd.length, 1);
  assert.equal(fwd[0].top, 9);
  assert.equal(fwd[0].bottom, 13);
  const rc = revcomp(seq);
  const rev = findSites(rc, false, ENZYME_BY_NAME.BsaI);
  assert.equal(rev.length, 1);
  assert.equal(rev[0].strand, -1);
  assert.equal(rev[0].top, rc.length - 13);
  assert.equal(rev[0].bottom, rc.length - 9);
});

test('digest: linear and circular fragment sizes add up', () => {
  const seq = 'A'.repeat(100) + 'GAATTC' + 'C'.repeat(200) + 'GGATCC' + 'T'.repeat(300);
  const cuts = ['EcoRI', 'BamHI'].flatMap((n) => findSites(seq, false, ENZYME_BY_NAME[n]));
  const lin = digest(seq, false, cuts);
  assert.equal(lin.length, 3);
  assert.equal(lin.reduce((n, f) => n + f.len, 0), seq.length);
  const circ = digest(seq, true, ['EcoRI', 'BamHI'].flatMap((n) => findSites(seq, true, ENZYME_BY_NAME[n])));
  assert.equal(circ.length, 2);
  assert.equal(circ.reduce((n, f) => n + f.len, 0), seq.length);
  const uncut = digest(seq, true, []);
  assert.equal(uncut.length, 1);
});

test('gel migration decreases with size', () => {
  assert.ok(migration(500) > migration(5000));
  assert.ok(migration(100, 2) > migration(1000, 2));
});

test('editing keeps features on the same bases', () => {
  const doc = newDocument({ seq: 'AAAAACCCCCGGGGGTTTTT', circular: false });
  doc.features = [feat('C block', 5, 10), feat('T block', 15, 20)];
  replaceRange(doc, 10, 15, 'XX'); // replace the G block
  assert.equal(doc.seq, 'AAAAACCCCCXXTTTTT');
  assert.equal(subseq(doc.seq, doc.features[0].segments[0].start, doc.features[0].segments[0].end), 'CCCCC');
  assert.equal(subseq(doc.seq, doc.features[1].segments[0].start, doc.features[1].segments[0].end), 'TTTTT');
  replaceRange(doc, 0, 0, 'GG'); // insertion at the start shifts everything
  assert.equal(doc.features[0].segments[0].start, 7);
  replaceRange(doc, 6, 9, ''); // delete into the start of the C block
  assert.equal(featureSequence(doc, doc.features[0]), 'CCC');
});

test('rotating the origin and reverse-complementing keep features intact', () => {
  const doc = newDocument({ seq: 'ATGAAACCCGGGTTTTAG' + 'A'.repeat(12), circular: true });
  doc.features = [feat('orf', 0, 18, 1)];
  const before = featureSequence(doc, doc.features[0]);
  rotateOrigin(doc, 10);
  assert.equal(featureSequence(doc, doc.features[0]), before);
  assert.ok(doc.features[0].segments[0].end <= doc.features[0].segments[0].start, 'feature now wraps the origin');
  reverseComplementDoc(doc);
  assert.equal(featureSequence(doc, doc.features[0]), before);
  assert.equal(doc.features[0].strand, -1);
});

test('ORF finder finds a forward and a reverse ORF', () => {
  const orf = 'ATG' + 'GCT'.repeat(110) + 'TAA';
  const seq = 'CC' + orf + 'CCCC' + revcomp(orf) + 'CC';
  const orfs = findOrfs(seq, false, 100);
  assert.equal(orfs.filter((o) => o.strand === 1 && o.codons === 111).length, 1);
  assert.equal(orfs.filter((o) => o.strand === -1 && o.codons === 111).length, 1);
  const r = orfs.find((o) => o.strand === -1);
  assert.equal(revcomp(seq.slice(r.start, r.end)), orf);
});

test('primer binding sites, design and PCR', () => {
  const left = 'GATTACAGATTACACCGGTTAAGGCC';
  const middle = 'T'.repeat(30) + 'GAATTC' + 'A'.repeat(30);
  const right = 'CCTTGGAACCTTGGCATCGATCGTAG';
  // flanks are T so the primer's 5' tail (…CC) can't happen to anneal
  const template = newDocument({ seq: 'T'.repeat(40) + left + middle + right + 'T'.repeat(40), circular: false });
  template.features = [feat('site', 40 + left.length + 30, 40 + left.length + 36)];
  const fwd = 'GGATCC' + left; // with a BamHI tail
  const rev = revcomp(right);
  const sites = findBindingSites(template.seq, false, fwd);
  assert.equal(sites.length, 1);
  assert.equal(sites[0].strand, 1);
  assert.equal(sites[0].tail, 'GGATCC');
  const { doc } = simulatePCR(template, fwd, rev);
  assert.equal(doc.seq, 'GGATCC' + left + middle + right);
  assert.equal(doc.features.length, 3);
  assert.equal(doc.features.filter((f) => f.type === 'primer_bind').length, 2);
  assert.equal(featureSequence(doc, doc.features.find((f) => f.type === 'misc_feature')), 'GAATTC');
  const design = designPrimers(template.seq, false, 40, 40 + left.length + middle.length + right.length, { targetTm: 58 });
  assert.ok(template.seq.includes(design.forward.seq));
  assert.ok(template.seq.includes(revcomp(design.reverse.seq)));
});

test('PCR across the origin of a circular template', () => {
  const seq = 'ACGTTGCAAGGCTTAACCGG' + 'A'.repeat(200) + 'TTGGCCAAGGTTCCAAGGTT';
  const doc = newDocument({ seq, circular: true });
  const fwd = seq.slice(-20, -2); // binds near the end
  const rev = revcomp(seq.slice(2, 20)); // binds near the start
  const { doc: product } = simulatePCR(doc, fwd, rev);
  assert.equal(product.seq, seq.slice(-20, -2) + seq.slice(-2) + seq.slice(0, 20));
});

test('XML round trip keeps attributes with markup in them', () => {
  const x = '<?xml version="1.0"?><Notes><Description>&lt;html&gt;hi&lt;/html&gt;</Description><R a="x &amp; y" b="&lt;b>"/></Notes>';
  const root = parseXML(x);
  assert.equal(root.children[0].text, '<html>hi</html>');
  assert.equal(root.children[1].attrs.b, '<b>');
  const again = parseXML(toXML(root));
  assert.deepEqual(again, root);
  assert.equal(toXML(el('A', { n: '1' })), '<A n="1"/>');
});

test('GenBank locations', () => {
  assert.deepEqual(parseLocation('complement(join(10..20,30..40))', 100), { strand: -1, segments: [{ start: 9, end: 20 }, { start: 29, end: 40 }] });
  assert.deepEqual(parseLocation('join(90..100,1..5)', 100), { strand: 1, segments: [{ start: 89, end: 5 }] });
  assert.deepEqual(parseLocation('<1..>50', 100), { strand: 1, segments: [{ start: 0, end: 50, partialStart: true, partialEnd: true }] });
  const f = { strand: -1, segments: [{ start: 89, end: 5 }] };
  assert.equal(formatLocation(f, 100), 'complement(join(90..100,1..5))');
});

test('GenBank and FASTA round trips', () => {
  const doc = newDocument({ name: 'pTest', seq: 'ATGC'.repeat(40), circular: true });
  doc.features = [feat('thing', 4, 20, -1), feat('wrap', 150, 10, 1)];
  doc.notes.description = 'A test plasmid';
  const [back] = readGenBank(writeGenBank(doc));
  assert.equal(back.seq, doc.seq);
  assert.equal(back.circular, true);
  assert.equal(back.features.length, 2);
  assert.deepEqual(back.features.map((f) => [f.strand, f.segments.map((s) => [s.start, s.end])]), [[-1, [[4, 20]]], [1, [[150, 10]]]]);
  const [fa] = readFasta(writeFasta(doc));
  assert.equal(fa.seq, doc.seq);
  assert.equal(fa.circular, true);
});

test('featuresInRegion handles wrap-around', () => {
  const doc = newDocument({ seq: 'A'.repeat(100), circular: true });
  doc.features = [feat('wrap', 95, 5)];
  const out = featuresInRegion(doc, 90, 20, 3);
  assert.equal(out.length, 1);
  assert.deepEqual([out[0].segments[0].start, out[0].segments[0].end], [8, 18]);
  assert.equal(rangeLength({ start: 95, end: 5 }, 100), 10);
});
