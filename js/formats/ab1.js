// ABI trace (.ab1) reader.
//
// The container is a small directory of 28-byte tag entries followed by data
// blocks. Data of four bytes or fewer is stored in the entry itself, in place
// of the offset. Directory offsets differ between writers, so the reader first
// tries the standard header offsets and then scans for a valid run of entries.
// Only the tags needed for a chromatogram are interpreted.
//
// Element types (ABIF): 1 byte, 2 char, 3 word (u16), 4 short (i16),
// 5 long (i32), 7 float, 8 double; 18 pString, 19 cString, 1023 and up are
// directories and user types, kept as opaque bytes.

const DEFAULT_SIZE = { 1: 1, 2: 1, 3: 2, 4: 2, 5: 4, 6: 8, 7: 4, 8: 8, 10: 4, 11: 4, 12: 10, 13: 1, 18: 1, 19: 1 };

function tagName(buf, p) {
  return String.fromCharCode(buf[p], buf[p + 1], buf[p + 2], buf[p + 3]);
}

function plainChars(buf, start, count) {
  let out = '';
  for (let i = 0; i < count; i++) {
    const c = buf[start + i];
    if (c === 0) break;
    if (c >= 32 && c < 127) out += String.fromCharCode(c);
  }
  return out;
}

function readEntry(dv, p, length) {
  if (p + 28 > length) return null;
  const view = new Uint8Array(dv.buffer, dv.byteOffset, dv.byteLength);
  const name = tagName(view, p);
  if (!/^[\x21-\x7e]{4}$/.test(name)) return null; // tag names are four printable characters, such as 'Dye#'
  const number = dv.getUint32(p + 4, false);
  const type = dv.getUint16(p + 8, false);
  const elemSize = dv.getUint16(p + 10, false);
  const count = dv.getUint32(p + 12, false);
  const dataSize = dv.getUint32(p + 16, false);
  // small values live in the entry, where the offset would be
  const dataOffset = dataSize <= 4 ? p + 20 : dv.getUint32(p + 20, false);
  const size = elemSize || DEFAULT_SIZE[type] || 1;
  if (dataOffset > length || dataOffset + dataSize > length) return null;
  return { name, number, key: `${name}${number}`, type, elemSize: size, count, dataSize, dataOffset };
}

function parseDirectory(dv, start, length) {
  const entries = [];
  for (let p = start; p + 28 <= length; p += 28) {
    const e = readEntry(dv, p, length);
    if (!e) break;
    entries.push(e);
    if (entries.length > 512) break;
  }
  return entries;
}

function findDirectory(buf, dv, length) {
  const candidates = [];
  // standard-ish locations for the directory offset/count in the 28-byte header
  for (const off of [10, 14, 18, 22, 26, 6]) {
    if (off + 4 > length) continue;
    const start = dv.getUint32(off, false);
    if (start < 28 || start + 28 > length) continue;
    const entries = parseDirectory(dv, start, length);
    if (entries.length >= 2) candidates.push(entries);
  }
  // fallback: scan for the first plausible entry run
  for (let p = 28; p + 28 * 2 <= length; p += 4) {
    const first = readEntry(dv, p, length);
    if (!first) continue;
    const entries = parseDirectory(dv, p, length);
    if (entries.length >= 2) { candidates.push(entries); break; }
  }
  candidates.sort((a, b) => b.length - a.length);
  return candidates[0] ?? [];
}

function values(buf, dv, entry) {
  const out = [];
  const base = entry.dataOffset;
  const n = Math.min(entry.count, Math.floor(entry.dataSize / Math.max(1, entry.elemSize)));
  for (let i = 0; i < n; i++) {
    const p = base + i * entry.elemSize;
    if (entry.type === 1 || entry.type === 2 || entry.type === 13) out.push(buf[p]);
    else if (entry.type === 3) out.push(dv.getUint16(p, false));
    else if (entry.type === 4) out.push(dv.getInt16(p, false));
    else if (entry.type === 5) out.push(dv.getInt32(p, false));
    else if (entry.type === 7) out.push(dv.getFloat32(p, false));
    else if (entry.type === 6 || entry.type === 8) out.push(dv.getFloat64(p, false));
  }
  return out;
}

function pleasant(buf, entry) {
  if (!entry) return '';
  if (entry.type === 2 || entry.type === 1) {
    return plainChars(buf, entry.dataOffset, entry.dataSize || entry.count);
  }
  // pString: first byte is the length
  const len = buf[entry.dataOffset];
  if (len != null && len + 1 <= entry.dataSize) return plainChars(buf, entry.dataOffset + 1, len);
  return plainChars(buf, entry.dataOffset, entry.dataSize);
}

export function readAB1(buffer, fileName = 'trace.ab1') {
  const buf = buffer instanceof Uint8Array ? buffer : new Uint8Array(buffer);
  if (buf.length < 32) throw new Error('Trace file is too small to be an ABI trace');
  const magic = String.fromCharCode(buf[0], buf[1], buf[2], buf[3]);
  if (magic !== 'ABIF') throw new Error('Not an ABI trace (missing ABIF header)');
  const dv = new DataView(buf.buffer, buf.byteOffset, buf.byteLength);
  const entries = findDirectory(buf, dv, buf.length);
  if (!entries.length) throw new Error('ABI trace has no readable tag directory');
  const byKey = new Map(entries.map((e) => [e.key, e]));
  const pick = (...names) => {
    for (const name of names) {
      for (const [key, e] of byKey) if (key.startsWith(name)) return e;
    }
    return null;
  };

  const basesEntry = pick('PBAS2', 'PBAS1', 'FWO_1', 'BASE');
  const qualEntry = pick('PCON2', 'PCON1');
  const locEntry = pick('PLOC2', 'PLOC1');
  const nameEntry = pick('SMPL1', 'SMPL');

  // The four dye channels are separate DATA tags (9-12 are the analysed ones,
  // 1-4 the raw), in the order FWO_1 gives: 'GATC' on most instruments.
  const chan = { A: new Uint16Array(0), C: new Uint16Array(0), G: new Uint16Array(0), T: new Uint16Array(0) };
  const order = (pleasant(buf, byKey.get('FWO_1')).toUpperCase().replace(/[^ACGT]/g, '') + 'ACGT').slice(0, 4);
  const tags = [9, 10, 11, 12].every((n) => byKey.has(`DATA${n}`)) ? [9, 10, 11, 12]
    : [1, 2, 3, 4].every((n) => byKey.has(`DATA${n}`)) ? [1, 2, 3, 4] : null;
  let sampleCount = 0;
  if (tags) {
    tags.forEach((n, c) => {
      const raw = values(buf, dv, byKey.get(`DATA${n}`));
      chan[order[c]] = Uint16Array.from(raw, (v) => Math.max(0, v)); // a trace can dip a little below zero
      sampleCount = Math.max(sampleCount, raw.length);
    });
  } else {
    // one tag carrying all four channels end to end
    const data = pick('DATA9') ?? pick('DATA');
    if (data) {
      const raw = values(buf, dv, data);
      sampleCount = Math.floor(raw.length / 4);
      ['A', 'C', 'G', 'T'].forEach((base, c) => {
        const arr = new Uint16Array(sampleCount);
        for (let i = 0; i < sampleCount; i++) arr[i] = raw[c * sampleCount + i] ?? 0;
        chan[base] = arr;
      });
    }
  }
  const bases = pleasant(buf, basesEntry).replace(/[^ACGTN]/gi, '').toUpperCase();
  const confidence = qualEntry ? values(buf, dv, qualEntry) : null;
  // peak positions are stored as signed shorts but are sample numbers, so read them as unsigned
  const positions = locEntry ? values(buf, dv, locEntry).map((v) => v & 0xffff) : [];
  if (positions.length && positions.length === bases.length && positions[0] > 0) {
    // ABI PLOC is 1-based in practice; a first peak at 0 is already 0-based.
    for (let i = 0; i < positions.length; i++) positions[i] -= 1;
  }
  const name = nameEntry ? pleasant(buf, nameEntry).replace(/\0+$/, '') : fileName.replace(/\.(ab1|abi)$/i, '');
  if (!sampleCount && !bases.length) throw new Error('ABI trace has no samples or bases');
  return {
    format: 'AB1',
    name,
    bases,
    positions,
    confidence,
    qualities: confidence,
    samples: chan,
    sampleCount: sampleCount || chan.A.length,
    tags: entries.map((e) => e.key),
  };
}
