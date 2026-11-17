// BAM reader, for the aligned reads stored inside a .dna file.
//
// BAM is a BGZF stream: a run of gzip members, each carrying its own compressed
// size in an extra field. The blocks are split on that size and inflated one at
// a time, rather than relying on the platform to decompress concatenated gzip
// members. Works in the browser (DecompressionStream) and in Node (zlib).

import { parseXML, children } from '../xml.js';

const SEQ_CHARS = '=ACMGRSVTWYHKDBN';
const CIGAR_OPS = 'MIDNSHP=X';

async function inflateRaw(bytes) {
  if (typeof DecompressionStream === 'function') {
    const stream = new Blob([bytes]).stream().pipeThrough(new DecompressionStream('deflate-raw'));
    return new Uint8Array(await new Response(stream).arrayBuffer());
  }
  const { inflateRawSync } = await import('node:zlib');
  return new Uint8Array(inflateRawSync(bytes));
}

// Split a BGZF stream into its deflate payloads and inflate each one.
export async function inflateBgzf(bytes) {
  const parts = [];
  let total = 0;
  let i = 0;
  while (i < bytes.length) {
    if (bytes[i] !== 0x1f || bytes[i + 1] !== 0x8b) throw new Error(`Not a gzip block at byte ${i}`);
    const flg = bytes[i + 3];
    if (!(flg & 4)) throw new Error('Not BGZF: the gzip block has no extra field');
    const xlen = bytes[i + 10] | (bytes[i + 11] << 8);
    let bsize = null;
    let p = i + 12;
    const end = p + xlen;
    while (p < end) {
      const si1 = bytes[p];
      const si2 = bytes[p + 1];
      const slen = bytes[p + 2] | (bytes[p + 3] << 8);
      if (si1 === 66 && si2 === 67 && slen === 2) bsize = (bytes[p + 4] | (bytes[p + 5] << 8)) + 1;
      p += 4 + slen;
    }
    if (!bsize) throw new Error('BGZF block without a size field');
    const payload = bytes.subarray(i + 12 + xlen, i + bsize - 8);
    if (payload.length) {
      const out = await inflateRaw(payload);
      parts.push(out);
      total += out.length;
    }
    i += bsize;
  }
  const joined = new Uint8Array(total);
  let o = 0;
  for (const part of parts) { joined.set(part, o); o += part.length; }
  return joined;
}

class Cursor {
  constructor(bytes) {
    this.b = bytes;
    this.v = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
    this.i = 0;
  }

  i32() { const n = this.v.getInt32(this.i, true); this.i += 4; return n; }
  u32() { const n = this.v.getUint32(this.i, true); this.i += 4; return n; }
  u16() { const n = this.v.getUint16(this.i, true); this.i += 2; return n; }
  u8() { const n = this.b[this.i]; this.i += 1; return n; }
  str(n) { const s = new TextDecoder().decode(this.b.subarray(this.i, this.i + n)); this.i += n; return s; }
}

function decodeSeq(cur, length) {
  let s = '';
  for (let k = 0; k < Math.ceil(length / 2); k++) {
    const byte = cur.u8();
    s += SEQ_CHARS[byte >> 4];
    if (s.length < length) s += SEQ_CHARS[byte & 15];
  }
  return s;
}

function decodeCigar(cur, n) {
  const ops = [];
  for (let k = 0; k < n; k++) {
    const raw = cur.u32();
    ops.push([raw >>> 4, CIGAR_OPS[raw & 15]]);
  }
  return ops;
}

// How far along the reference a read reaches.
export function referenceSpan(cigar) {
  let n = 0;
  for (const [len, op] of cigar) if ('MDN=X'.includes(op)) n += len;
  return n;
}

export async function readBam(bytes) {
  const raw = await inflateBgzf(bytes instanceof Uint8Array ? bytes : new Uint8Array(bytes));
  const cur = new Cursor(raw);
  if (cur.str(4) !== 'BAM\u0001') throw new Error('Not a BAM file');
  const headerText = cur.str(cur.i32());
  const refs = [];
  const nRef = cur.i32();
  for (let k = 0; k < nRef; k++) {
    const name = cur.str(cur.i32()).replace(/\0+$/, '');
    refs.push({ name, length: cur.i32() });
  }
  const reads = [];
  while (cur.i + 4 <= raw.length) {
    const blockSize = cur.i32();
    if (blockSize <= 0 || cur.i + blockSize > raw.length) break;
    const next = cur.i + blockSize;
    const refId = cur.i32();
    const pos = cur.i32();
    const lReadName = cur.u8();
    const mapq = cur.u8();
    cur.u16(); // bin, which indexes a file we are not seeking in
    const nCigar = cur.u16();
    const flag = cur.u16();
    const lSeq = cur.i32();
    cur.i32(); cur.i32(); cur.i32(); // mate reference, mate position, template length
    const name = cur.str(lReadName).replace(/\0+$/, '');
    const cigar = decodeCigar(cur, nCigar);
    const seq = decodeSeq(cur, lSeq);
    const qual = raw.subarray(cur.i, cur.i + lSeq);
    cur.i = next; // skip the optional tags
    reads.push({
      name, refId, pos, mapq, flag, cigar, seq,
      qual: qual[0] === 0xff ? null : Array.from(qual),
      reverse: !!(flag & 16),
      unmapped: !!(flag & 4),
      refLength: referenceSpan(cigar),
      reference: refs[refId]?.name ?? null,
    });
  }
  return { headerText, refs, reads };
}

// Pull the aligned reads out of an opened .dna document. The reads live in the
// alignment segments; their display names and trim ranges come from the
// alignable-sequences list beside them, keyed by the read name.
export async function readAlignments(doc) {
  const segs = doc.source?.segments ?? [];
  const meta = new Map();
  for (const s of segs) {
    if (s.type !== 0x11) continue;
    try {
      for (const el of children(parseXML(new TextDecoder().decode(s.data)), 'Sequence')) {
        meta.set(el.attrs.ID, el.attrs);
      }
    } catch { /* an unreadable list just means unnamed reads */ }
  }
  const reads = [];
  for (const s of segs) {
    if (s.type !== 0x1b) continue;
    try {
      const bam = await readBam(s.data);
      for (const r of bam.reads) {
        if (r.unmapped) continue;
        const m = meta.get(r.name);
        reads.push({ ...r, label: m?.name || r.name, trimmedRange: m?.trimmedRange ?? null, ...alignToReference(r) });
      }
    } catch (e) {
      /** @type {any} */ (reads).errors = [...(/** @type {any} */ (reads).errors ?? []), e.message];
    }
  }
  reads.sort((a, b) => a.pos - b.pos || a.label.localeCompare(b.label));
  return reads.length ? reads : null;
}

// One character per reference position the read covers: the read's base, or
// '-' where the read is missing bases the reference has. Insertions are noted
// against the position they follow.
export function alignToReference(read) {
  const chars = [];
  const insertions = new Map();
  let q = 0;
  for (const [len, op] of read.cigar) {
    if (op === 'M' || op === '=' || op === 'X') {
      for (let k = 0; k < len; k++) chars.push(read.seq[q + k] ?? 'N');
      q += len;
    } else if (op === 'D' || op === 'N') {
      for (let k = 0; k < len; k++) chars.push('-');
    } else if (op === 'I') {
      insertions.set(read.pos + chars.length, read.seq.slice(q, q + len));
      q += len;
    } else if (op === 'S') {
      q += len; // soft-clipped: present in the read, not aligned
    }
    // 'H' and 'P' consume neither the read nor the reference
  }
  return { start: read.pos, chars, insertions };
}
