// ZTR reader, for the sequencing traces stored inside a .dna file.
//
// A ZTR file is a header and a series of chunks. Each chunk's data is wrapped
// in a chain of transforms, innermost last, each announcing itself with a
// leading format byte; the chain is unwound until the raw data is reached.
// The transforms seen in practice are zlib, run-length, "follow" (predict each
// byte from the one before it), 16/32-bit-to-8-bit packing, and delta.
//
// The raw data keeps its leading format byte, padded out to the width of its
// values, so that A, C, G and T samples stay two-byte aligned. Each chunk
// skips that leading value itself.

const MAGIC = [0xae, 0x5a, 0x54, 0x52, 0x0d, 0x0a, 0x1a, 0x0a];
const RAW = 0, RLE = 1, ZLIB = 2;
const DELTA1 = 64, DELTA2 = 65, DELTA4 = 66;
const TO8_16 = 70, TO8_32 = 71, FOLLOW1 = 72;

async function inflateZlib(bytes) {
  if (typeof DecompressionStream === 'function') {
    const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('deflate'));
    return new Uint8Array(await new Response(stream).arrayBuffer());
  }
  const { inflateSync } = await import('node:zlib');
  return new Uint8Array(inflateSync(bytes));
}

function unRle(d) {
  const guard = d[5];
  const out = [];
  for (let i = 6; i < d.length;) {
    const b = d[i++];
    if (b !== guard) { out.push(b); continue; }
    const count = d[i++];
    if (count === 0) { out.push(guard); continue; }
    const value = d[i++];
    for (let k = 0; k < count; k++) out.push(value);
  }
  return Uint8Array.from(out);
}

// Each byte is stored as how far it falls short of the byte that usually
// follows the one before it. The first byte is stored as it is.
function unFollow(d) {
  const table = d.subarray(1, 257);
  const body = d.subarray(257);
  const out = new Uint8Array(body.length);
  if (!body.length) return out;
  out[0] = body[0];
  for (let i = 1; i < body.length; i++) out[i] = (table[out[i - 1]] - body[i]) & 0xff;
  return out;
}

// Values that fit in a signed byte take one byte; the rest are escaped with
// -128 and written out in full, big-endian.
function unPack(d, size) {
  const out = [];
  for (let i = 1; i < d.length;) {
    const b = d[i++];
    if (b === 0x80) {
      for (let k = 0; k < size && i + k < d.length; k++) out.push(d[i + k]);
      i += size;
    } else {
      const v = b > 127 ? b - 256 : b;
      for (let k = size - 1; k >= 0; k--) out.push((v >> (8 * k)) & 0xff);
    }
  }
  return Uint8Array.from(out);
}

// Undo repeated differencing, wrapping within the value width exactly as the
// writer did.
function unDelta(d, size) {
  const level = d[1];
  const header = size === 4 ? 4 : 2;
  const body = d.subarray(header);
  const n = Math.floor(body.length / size);
  const read = new DataView(body.buffer, body.byteOffset, body.byteLength);
  const vals = new Array(n);
  for (let k = 0; k < n; k++) {
    vals[k] = size === 1 ? body[k] : size === 2 ? read.getUint16(k * 2, false) : read.getUint32(k * 4, false);
  }
  const mask = size === 1 ? 0xff : size === 2 ? 0xffff : 0xffffffff;
  for (let pass = 0; pass < level; pass++) {
    let run = 0;
    for (let k = 0; k < n; k++) {
      run = size === 4 ? ((run + vals[k]) >>> 0) : ((run + vals[k]) & mask);
      vals[k] = run;
    }
  }
  const out = new Uint8Array(n * size);
  const write = new DataView(out.buffer);
  for (let k = 0; k < n; k++) {
    if (size === 1) out[k] = vals[k];
    else if (size === 2) write.setUint16(k * 2, vals[k], false);
    else write.setUint32(k * 4, vals[k], false);
  }
  return out;
}

async function decodeChain(data) {
  let d = data;
  for (let guard = 0; guard < 12; guard++) {
    const form = d[0];
    if (form === RAW) return d;              // leading marker kept, see the note above
    if (form === ZLIB) d = await inflateZlib(d.subarray(5));
    else if (form === RLE) d = unRle(d);
    else if (form === FOLLOW1) d = unFollow(d);
    else if (form === TO8_16) d = unPack(d, 2);
    else if (form === TO8_32) d = unPack(d, 4);
    else if (form === DELTA1) d = unDelta(d, 1);
    else if (form === DELTA2) d = unDelta(d, 2);
    else if (form === DELTA4) d = unDelta(d, 4);
    else throw new Error(`Unsupported ZTR encoding 0x${form.toString(16)}`);
  }
  throw new Error('ZTR encodings nested too deeply');
}

export async function readZtr(bytes) {
  const d = bytes instanceof Uint8Array ? bytes : new Uint8Array(bytes);
  for (let i = 0; i < MAGIC.length; i++) if (d[i] !== MAGIC[i]) throw new Error('Not a ZTR trace');
  const head = new DataView(d.buffer, d.byteOffset, d.byteLength);
  const trace = { samples: null, bases: '', positions: [], confidence: null, text: {} };

  for (let i = 10; i + 8 <= d.length;) {
    const type = String.fromCharCode(d[i], d[i + 1], d[i + 2], d[i + 3]);
    const metaLen = head.getUint32(i + 4, false);
    const p = i + 8 + metaLen;
    if (p + 4 > d.length) break;
    const dataLen = head.getUint32(p, false);
    const body = d.subarray(p + 4, p + 4 + dataLen);
    i = p + 4 + dataLen;

    let raw;
    try { raw = await decodeChain(body); } catch { continue; }
    const v = new DataView(raw.buffer, raw.byteOffset, raw.byteLength);

    if (type === 'SMP4') {
      const n = Math.floor((Math.floor(raw.length / 2) - 1) / 4);
      const chan = {};
      ['A', 'C', 'G', 'T'].forEach((name, k) => {
        const arr = new Uint16Array(n);
        for (let s = 0; s < n; s++) arr[s] = v.getUint16(2 + (k * n + s) * 2, false);
        chan[name] = arr;
      });
      trace.samples = chan;
    } else if (type === 'BASE') {
      trace.bases = new TextDecoder().decode(raw.subarray(1)).replace(/\0+$/, '');
    } else if (type === 'BPOS') {
      const n = Math.floor(raw.length / 4);
      for (let k = 1; k < n; k++) trace.positions.push(v.getUint32(k * 4, false));
    } else if (type === 'CNF4' || type === 'CNF1') {
      trace.confidence = Array.from(raw.subarray(1));
    } else if (type === 'TEXT') {
      const parts = new TextDecoder().decode(raw.subarray(1)).split('\0');
      for (let k = 0; k + 1 < parts.length; k += 2) if (parts[k]) trace.text[parts[k]] = parts[k + 1];
    }
  }

  if (!trace.samples) throw new Error('The trace has no sample data');
  trace.sampleCount = trace.samples.A.length;
  trace.name = trace.text.NAME ?? '';
  return trace;
}

// Pull the traces out of an opened .dna document. Each sits in its own segment
// behind a 4-byte id and a nested chunk header.
export async function readTraces(doc) {
  const out = [];
  for (const seg of doc.source?.segments ?? []) {
    if (seg.type !== 0x10) continue;
    const d = seg.data;
    const v = new DataView(d.buffer, d.byteOffset, d.byteLength);
    const id = String(v.getUint32(0, false));
    const innerLen = v.getUint32(5, false);
    try {
      out.push({ id, ...(await readZtr(d.subarray(9, 9 + innerLen))) });
    } catch (e) {
      out.push({ id, error: e.message });
    }
  }
  return out;
}
