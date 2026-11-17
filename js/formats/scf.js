// SCF trace reader (.scf), versions 2 and 3.
//
// Samples are four interleaved channels; base records carry a peak index and a
// probability, which becomes a Phred score when an explicit quality block is
// absent. Version 3's extended header supplies the sample width.

function u32(dv, p) {
  if (p + 4 > dv.byteLength) return 0;
  return dv.getUint32(p, false);
}

export function readSCF(buffer, fileName = 'trace.scf') {
  const buf = buffer instanceof Uint8Array ? buffer : new Uint8Array(buffer);
  if (buf.length < 56) throw new Error('Trace file is too small to be an SCF trace');
  if (String.fromCharCode(buf[0], buf[1], buf[2], buf[3]) !== '.scf') throw new Error('Not an SCF trace (missing .scf header)');
  const dv = new DataView(buf.buffer, buf.byteOffset, buf.byteLength);
  const version = u32(dv, 4);
  const samples = u32(dv, 8);
  const sampleOffset = u32(dv, 12);
  const baseCount = u32(dv, 16);
  const basesOffset = u32(dv, 20);
  const sampleSize = version >= 3 ? u32(dv, 32) : 1;
  if (!samples || sampleOffset + samples * 4 * sampleSize > buf.length) throw new Error('SCF sample block is truncated');
  if (!baseCount || basesOffset + baseCount * 4 > buf.length) throw new Error('SCF base block is truncated');

  const chan = { A: new Uint16Array(samples), C: new Uint16Array(samples), G: new Uint16Array(samples), T: new Uint16Array(samples) };
  const order = ['A', 'C', 'G', 'T'];
  for (let i = 0; i < samples; i++) {
    for (let c = 0; c < 4; c++) {
      const p = sampleOffset + (i * 4 + c) * sampleSize;
      const value = sampleSize === 2 ? dv.getUint16(p, false) : buf[p];
      chan[order[c]][i] = value;
    }
  }
  let bases = '';
  const positions = [];
  const confidence = [];
  for (let i = 0; i < baseCount; i++) {
    const p = basesOffset + i * 4;
    const peak = buf[p];
    const prob = buf[p + 1];
    const base = String.fromCharCode(buf[p + 2]).toUpperCase();
    positions.push(peak);
    const q = prob === 0 ? 99 : Math.min(99, Math.round(-10 * Math.log10(Math.max(1e-4, prob / 255))));
    confidence.push(q);
    bases += /[ACGTN]/.test(base) ? base : 'N';
  }
  return {
    format: 'SCF',
    name: fileName.replace(/\.scf$/i, ''),
    bases,
    positions,
    confidence,
    qualities: confidence,
    samples: chan,
    sampleCount: samples,
    version,
  };
}
