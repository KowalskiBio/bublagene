export function writePString(str, buf, offset) {
  buf[offset] = str.length;
  for (let i = 0; i < str.length; i++) buf[offset + 1 + i] = str.charCodeAt(i);
}
export function makeAB1({ bases = 'ACGTACGTAC', qualities = [5, 5, 5, 30, 30, 30, 30, 30, 5, 5], name = 'sample1' } = {}) {
  const n = bases.length;
  const tags = [];
  const dirOffset = 128;
  const entrySize = 28;
  let dataOffset = dirOffset + entrySize * 5;
  const add = (name, number, type, elemSize, count, dataSize) => {
    const tag = { name, number, type, elemSize, count, dataSize, dataOffset };
    tags.push(tag);
    dataOffset += dataSize;
  };
  const dataSize = n * 4 * 2;
  add('DATA', 9, 3, 2, n * 4, dataSize);
  add('PBAS', 2, 2, 1, n, n);
  add('PCON', 2, 1, 1, n, n);
  add('PLOC', 2, 3, 2, n, n * 2);
  add('SMPL', 1, 2, 1, name.length + 1, name.length + 1);
  const total = dataOffset;
  const buf = new Uint8Array(total);
  const dv = new DataView(buf.buffer);
  buf.set([0x41, 0x42, 0x49, 0x46], 0); // ABIF
  dv.setUint16(4, 101, false);
  dv.setUint32(10, dirOffset, false);
  dv.setUint32(14, tags.length, false);
  for (const [k, t] of tags.entries()) {
    const p = dirOffset + k * entrySize;
    for (let c = 0; c < 4; c++) buf[p + c] = t.name.charCodeAt(c);
    dv.setUint32(p + 4, t.number, false);
    dv.setUint16(p + 8, t.type, false);
    dv.setUint16(p + 10, t.elemSize, false);
    dv.setUint32(p + 12, t.count, false);
    dv.setUint32(p + 16, t.dataSize, false);
    dv.setUint32(p + 20, t.dataOffset, false);
    dv.setUint32(p + 24, 0, false);
  }
  const byKey = Object.fromEntries(tags.map((t) => [`${t.name}${t.number}`, t]));
  const data = byKey.DATA9;
  for (let i = 0; i < n; i++) {
    for (let c = 0; c < 4; c++) dv.setUint16(data.dataOffset + (c * n + i) * 2, 100 + c * 10 + i, false);
  }
  for (let i = 0; i < n; i++) {
    buf[byKey.PBAS2.dataOffset + i] = bases.charCodeAt(i);
    buf[byKey.PCON2.dataOffset + i] = qualities[i];
    dv.setUint16(byKey.PLOC2.dataOffset + i * 2, 10 + i * 5 + 1, false); // 1-based
  }
  writePString(name, buf, byKey.SMPL1.dataOffset);
  return buf;
}


// A trace laid out the way real instruments write it: a root directory entry in
// the header, four separate DATA tags (9-12) in the order FWO_1 gives, small
// values stored inside their directory entry, signed-short peak positions, and
// tags such as 'Dye#' whose names are not plain letters and digits.
export function makeStandardAB1({ bases = 'ACGTACGTAC', name = 'SEQ123_A1', order = 'GATC', spacing = 12 } = {}) {
  const n = bases.length;
  const samples = n * spacing + 20;
  const entries = [];
  const bytes = (b) => Uint8Array.from(b);
  const be16 = (arr) => { const o = new Uint8Array(arr.length * 2); const d = new DataView(o.buffer); arr.forEach((v, i) => d.setInt16(i * 2, v, false)); return o; };
  const add = (tag, number, type, elem, count, data) => entries.push({ tag, number, type, elem, count, data });
  add('AEPt', 1, 4, 2, 1, be16([2]));
  add('DATA', 9, 4, 2, samples, null); add('DATA', 10, 4, 2, samples, null);
  add('DATA', 11, 4, 2, samples, null); add('DATA', 12, 4, 2, samples, null);
  add('Dye#', 1, 4, 2, 1, be16([4]));
  add('FWO_', 1, 2, 1, 4, bytes([...order].map((c) => c.charCodeAt(0))));
  add('PBAS', 2, 2, 1, n, bytes([...bases].map((c) => c.charCodeAt(0))));
  add('PCON', 2, 2, 1, n, bytes(Array(n).fill(40)));
  add('PLOC', 2, 4, 2, n, be16(Array.from({ length: n }, (_, i) => 10 + i * spacing)));
  add('SMPL', 1, 18, 1, name.length + 1, bytes([name.length, ...[...name].map((c) => c.charCodeAt(0))]));
  // channel curves: a bump of height 1000 in the channel of the called base, centred on each peak
  const channels = { A: [], C: [], G: [], T: [] };
  for (const c of 'ACGT') channels[c] = new Array(samples).fill(0);
  for (let i = 0; i < n; i++) channels[bases[i]][10 + i * spacing] = 1000;
  for (const [k, e] of entries.entries()) if (e.tag === 'DATA') e.data = be16(channels[order[e.number - 9]]);
  let offset = 128;
  for (const e of entries) {
    e.size = e.data.length;
    if (e.size > 4) { e.offset = offset; offset += e.size + (e.size % 2); }
  }
  const dirOffset = offset;
  const buf = new Uint8Array(dirOffset + entries.length * 28);
  const dv = new DataView(buf.buffer);
  buf.set([0x41, 0x42, 0x49, 0x46], 0);
  dv.setUint16(4, 101, false);
  const putEntry = (p, tag, number, type, elem, count, size, off) => {
    for (let c = 0; c < 4; c++) buf[p + c] = tag.charCodeAt(c);
    dv.setUint32(p + 4, number, false); dv.setUint16(p + 8, type, false); dv.setUint16(p + 10, elem, false);
    dv.setUint32(p + 12, count, false); dv.setUint32(p + 16, size, false); dv.setUint32(p + 20, off, false);
  };
  putEntry(6, 'tdir', 1, 1023, 28, entries.length, entries.length * 28, dirOffset);
  entries.forEach((e, k) => {
    const p = dirOffset + k * 28;
    putEntry(p, e.tag, e.number, e.type, e.elem, e.count, e.size, 0);
    if (e.size > 4) { dv.setUint32(p + 20, e.offset, false); buf.set(e.data, e.offset); } else buf.set(e.data, p + 20);
  });
  return buf;
}
