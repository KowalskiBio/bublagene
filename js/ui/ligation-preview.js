// A picture of a restriction ligation: the two digested pieces with their real
// sticky ends, and the molecule they make.

import { h, s, fmt } from './dom.js';
import { pieceFromFragment, restrictionClone, describeEnd } from '../cloning.js';
import { circularMap, linearMap } from './map.js';

const W = 300;
const PAD = 22;      // room for an overhang to stick out
const OH = 9;        // how far an overhang is drawn

// One piece as two strands. A 5′ overhang on the left is the top strand running
// past the bottom one; on the right it is the bottom strand that runs past.
function strands(root, piece, y, x0, x1, color) {
  const l = piece.left.type;
  const r = piece.right.type;
  const topL = l === '5' ? x0 - OH : x0;
  const botL = l === '3' ? x0 - OH : x0;
  const topR = r === '3' ? x1 + OH : x1;
  const botR = r === '5' ? x1 + OH : x1;
  root.append(s('rect', { x: topL, y, width: Math.max(2, topR - topL), height: 7, rx: 2, fill: color }));
  root.append(s('rect', { x: botL, y: y + 8, width: Math.max(2, botR - botL), height: 7, rx: 2, fill: color, opacity: 0.55 }));
}

function endLabel(root, end, x, y, anchor) {
  const names = end.enzyme ? ` · ${end.enzyme}` : '';
  root.append(s('text', { x, y, 'font-size': 10, fill: 'var(--muted)', 'text-anchor': anchor }, `${describeEnd({ type: end.type, oh: end.oh })}${names}`));
}

export function ligationPreview({ vector, insert, circular = true, orientation = 'auto' }) {
  let vp;
  let ip;
  let product;
  try {
    vp = pieceFromFragment(vector.doc, vector.fragment);
    ip = pieceFromFragment(insert.doc, insert.fragment);
    product = restrictionClone({
      backbone: { doc: vector.doc, fragment: vector.fragment },
      insert: { doc: insert.doc, fragment: insert.fragment },
      circular, orientation,
    });
  } catch {
    return null; // the panel already explains why these ends do not ligate
  }
  const box = h('div', { style: { margin: '10px 0' } });
  const longest = Math.max(vp.doc.seq.length, ip.doc.seq.length);
  const width = (n) => Math.max(36, ((W - 2 * PAD) * n) / longest);
  const root = s('svg', { viewBox: `0 0 ${W} 112`, width: '100%', role: 'img', 'aria-label': 'Digested pieces before ligation' });
  let y = 8;
  for (const { piece, color, name } of [{ piece: vp, color: 'var(--muted)', name: 'vector' }, { piece: ip, color: '#e0884c', name: 'insert' }]) {
    const w = width(piece.doc.seq.length);
    const x0 = PAD;
    root.append(s('text', { x: PAD, y: y + 1, 'font-size': 10.5, fill: 'var(--ink)', 'font-weight': 600 }, `${name}: ${piece.doc.name} · ${fmt(piece.doc.seq.length)} bp`));
    strands(root, piece, y + 6, x0, x0 + w, color);
    // labels sit at the edges of the drawing, so a short piece cannot make them collide
    endLabel(root, piece.left, PAD - OH, y + 33, 'start');
    endLabel(root, piece.right, W - PAD + OH, y + 33, 'end');
    y += 52;
  }
  box.append(root);
  const arrow = h('div', { class: 'sub', style: { textAlign: 'center', margin: '2px 0 4px' } }, `↓ ligation → ${fmt(product.doc.seq.length)} bp, insert ${product.orientation}`);
  box.append(arrow);
  const map = product.doc.circular ? circularMap(product.doc, {}) : linearMap(product.doc, {});
  map.removeAttribute('height');
  map.setAttribute('width', '100%');
  map.style.maxWidth = product.doc.circular ? '300px' : '100%';
  map.style.display = 'block';
  map.style.margin = '0 auto';
  box.append(map);
  return box;
}
