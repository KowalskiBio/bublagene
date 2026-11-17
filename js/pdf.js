// A minimal PDF writer: one page, a title line and one JPEG image scaled to fit.
// Enough to export a map or a plot without a PDF library.

const enc = new TextEncoder();
const bytes = (s) => enc.encode(s);

// PDF strings escape \ ( ) and must stay ASCII for the standard Helvetica font.
const pdfText = (s) => String(s).replace(/[^\x20-\x7e]/g, '?').replace(/([\\()])/g, '\\$1');

/**
 * @param {Uint8Array} jpeg a baseline JPEG
 * @param {{width:number,height:number,title?:string,subtitle?:string}} options pixel size of the image and optional heading
 */
export function imagePdf(jpeg, { width, height, title = '', subtitle = '' }) {
  const landscape = width > height * 1.15;
  const [pw, ph] = landscape ? [842, 595] : [595, 842]; // A4 in points
  const margin = 36;
  const head = title || subtitle ? 44 : 0;
  const maxW = pw - 2 * margin;
  const maxH = ph - 2 * margin - head;
  const scale = Math.min(maxW / width, maxH / height);
  const w = width * scale;
  const h = height * scale;
  const x = (pw - w) / 2;
  const y = ph - margin - head - h;

  const lines = [];
  if (title) lines.push(`BT /F1 16 Tf ${margin} ${ph - margin - 14} Td (${pdfText(title)}) Tj ET`);
  if (subtitle) lines.push(`BT /F1 10 Tf ${margin} ${ph - margin - 30} Td (${pdfText(subtitle)}) Tj ET`);
  lines.push(`q ${w.toFixed(2)} 0 0 ${h.toFixed(2)} ${x.toFixed(2)} ${y.toFixed(2)} cm /Im0 Do Q`);
  const content = lines.join('\n');

  const parts = [];
  const offsets = [];
  let length = 0;
  const push = (chunk) => { parts.push(chunk); length += chunk.length; };
  const object = (n, body, stream = null) => {
    offsets[n] = length;
    push(bytes(`${n} 0 obj\n${body}\n`));
    if (stream) { push(bytes('stream\n')); push(stream); push(bytes('\nendstream\n')); }
    push(bytes('endobj\n'));
  };

  push(bytes('%PDF-1.4\n'));
  object(1, '<< /Type /Catalog /Pages 2 0 R >>');
  object(2, '<< /Type /Pages /Kids [3 0 R] /Count 1 >>');
  object(3, `<< /Type /Page /Parent 2 0 R /MediaBox [0 0 ${pw} ${ph}] /Resources << /Font << /F1 5 0 R >> /XObject << /Im0 6 0 R >> >> /Contents 4 0 R >>`);
  const contentBytes = bytes(content);
  object(4, `<< /Length ${contentBytes.length} >>`, contentBytes);
  object(5, '<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>');
  object(6, `<< /Type /XObject /Subtype /Image /Width ${width} /Height ${height} /ColorSpace /DeviceRGB /BitsPerComponent 8 /Filter /DCTDecode /Length ${jpeg.length} >>`, jpeg);
  const xref = length;
  let table = `xref\n0 7\n0000000000 65535 f \n`;
  for (let n = 1; n <= 6; n++) table += `${String(offsets[n]).padStart(10, '0')} 00000 n \n`;
  push(bytes(`${table}trailer\n<< /Size 7 /Root 1 0 R >>\nstartxref\n${xref}\n%%EOF\n`));

  const out = new Uint8Array(length);
  let at = 0;
  for (const p of parts) { out.set(p, at); at += p.length; }
  return out;
}
