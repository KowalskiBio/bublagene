// Minimal XML parser/serializer, enough for the XML segments in .dna files.
// Works the same in the browser and in Node (no DOMParser needed).
// Node shape: { tag, attrs: {name: value}, children: [node], text }

const ENT = { lt: '<', gt: '>', amp: '&', quot: '"', apos: "'" };

export function decodeEntities(s) {
  return s.replace(/&(#x[0-9a-fA-F]+|#\d+|[a-zA-Z]+);/g, (m, e) => {
    if (e[0] === '#') {
      const code = e[1] === 'x' || e[1] === 'X' ? parseInt(e.slice(2), 16) : parseInt(e.slice(1), 10);
      return Number.isFinite(code) ? String.fromCodePoint(code) : m;
    }
    return ENT[e] ?? m;
  });
}

function tagEnd(str, i) {
  // index of the '>' closing the tag that starts at i, skipping quoted values
  let q = null;
  for (let j = i + 1; j < str.length; j++) {
    const c = str[j];
    if (q) { if (c === q) q = null; } else if (c === '"' || c === "'") q = c; else if (c === '>') return j;
  }
  throw new Error('Unterminated XML tag');
}

export function parseXML(str) {
  const doc = { tag: '#document', attrs: {}, children: [], text: '' };
  const stack = [doc];
  let i = 0;
  while (i < str.length) {
    const lt = str.indexOf('<', i);
    const textEnd = lt < 0 ? str.length : lt;
    if (textEnd > i) stack[stack.length - 1].text += decodeEntities(str.slice(i, textEnd));
    if (lt < 0) break;
    if (str.startsWith('<?', lt)) { i = str.indexOf('?>', lt) + 2; continue; }
    if (str.startsWith('<!--', lt)) { i = str.indexOf('-->', lt) + 3; continue; }
    if (str.startsWith('<![CDATA[', lt)) {
      const j = str.indexOf(']]>', lt);
      stack[stack.length - 1].text += str.slice(lt + 9, j);
      i = j + 3;
      continue;
    }
    if (str.startsWith('<!', lt)) { i = tagEnd(str, lt) + 1; continue; }
    const gt = tagEnd(str, lt);
    if (str[lt + 1] === '/') {
      if (stack.length > 1) stack.pop();
      i = gt + 1;
      continue;
    }
    let inner = str.slice(lt + 1, gt);
    const selfClosing = inner.endsWith('/');
    if (selfClosing) inner = inner.slice(0, -1);
    const m = /^([^\s/>]+)/.exec(inner);
    const node = { tag: m[1], attrs: {}, children: [], text: '' };
    const re = /([^\s=]+)\s*=\s*("([^"]*)"|'([^']*)')/g;
    let a;
    const rest = inner.slice(m[1].length);
    while ((a = re.exec(rest))) node.attrs[a[1]] = decodeEntities(a[3] ?? a[4] ?? '');
    stack[stack.length - 1].children.push(node);
    if (!selfClosing) stack.push(node);
    i = gt + 1;
  }
  const root = doc.children[0];
  if (!root) throw new Error('Empty XML');
  return root;
}

export function child(node, tag) {
  return node?.children.find((c) => c.tag === tag) ?? null;
}

export function children(node, tag) {
  return node ? node.children.filter((c) => c.tag === tag) : [];
}

export function el(tag, attrs = {}, kids = [], text = '') {
  const clean = {};
  for (const [k, v] of Object.entries(attrs)) if (v !== undefined && v !== null) clean[k] = String(v);
  return { tag, attrs: clean, children: kids, text };
}

const escText = (s) => s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const escAttr = (s) => escText(s).replace(/"/g, '&quot;');

export function toXML(node, { declaration = false } = {}) {
  const out = [];
  const walk = (n) => {
    const attrs = Object.entries(n.attrs).map(([k, v]) => ` ${k}="${escAttr(v)}"`).join('');
    const text = n.text ? escText(n.text) : '';
    if (!n.children.length && !text) { out.push(`<${n.tag}${attrs}/>`); return; }
    out.push(`<${n.tag}${attrs}>`);
    if (text) out.push(text);
    for (const c of n.children) walk(c);
    out.push(`</${n.tag}>`);
  };
  walk(node);
  return (declaration ? '<?xml version="1.0"?>' : '') + out.join('');
}

// Rich text is stored as "<html><body>…</body></html>"; strip to plain text.
export function htmlToText(s) {
  if (!s) return '';
  return decodeEntities(s.replace(/<br\s*\/?>/gi, '\n').replace(/<[^>]+>/g, '')).trim();
}
