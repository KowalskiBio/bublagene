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

