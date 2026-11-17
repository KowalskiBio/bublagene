// How a construct was made, as a tree: the product at the top, the pieces it was
// joined from beneath it, and, when those are open and were cloned themselves,
// what they were made from in turn.

import { h, fmt } from './dom.js';

const OPERATION = {
  'restriction-ligation': 'Restriction cloning', 'gibson-assembly': 'Gibson assembly', 'in-fusion-assembly': 'In-Fusion assembly',
  'golden-gate-assembly': 'Golden Gate assembly', 'site-specific-recombination': 'Recombination',
  'ta-topo-cloning': 'TA cloning', 'blunt-topo-cloning': 'Blunt cloning', 'directional-topo-cloning': 'Directional TOPO cloning',
};
const label = (op) => OPERATION[op] ?? String(op ?? 'cloning').replace(/-/g, ' ');

// A node for `doc`, with children for each input that is itself an open, cloned document.
export function buildTree(doc, docs, seen = new Set(), depth = 0) {
  const c = doc.cloning;
  const node = { id: doc.id, name: doc.name, size: doc.seq.length, circular: doc.circular, operation: c ? label(c.operation) : null, children: [], inputs: [] };
  if (!c || depth > 6 || seen.has(doc.id)) return node;
  seen.add(doc.id);
  c.inputs.forEach((i, k) => {
    const source = docs.find((d) => d.id !== doc.id && d.name === i.name && !seen.has(d.id));
    const junctions = c.junctions.filter((j) => j.leftInput === k || j.rightInput === k);
    const child = source ? buildTree(source, docs, seen, depth + 1) : { name: i.name, size: i.seqLen, children: [], inputs: [] };
    node.children.push({ ...child, role: i.role, orientation: i.orientation, enzymes: [...new Set(junctions.flatMap((j) => j.enzymes ?? []).filter(Boolean))], openId: source?.id ?? null });
  });
  return node;
}

function renderNode(node, activeId, onOpen) {
  const meta = [
    node.role, node.orientation && node.orientation !== 'forward' ? node.orientation : null,
    node.size != null ? `${fmt(node.size)} bp` : null, node.enzymes?.length ? node.enzymes.join('/') : null,
  ].filter(Boolean).join(' · ');
  const target = node.openId ?? node.id;
  const row = h('div', {
    class: 'treenode', 'aria-current': String(!!node.id && node.id === activeId),
    ...(target ? { role: 'button', tabindex: '0', title: 'Show this file', onclick: () => onOpen?.(target), onkeydown: (e) => { if (e.key === 'Enter') onOpen?.(target); } } : {}),
  },
  h('span', { class: 'treename' }, node.name),
  node.operation ? h('span', { class: 'treeop' }, node.operation) : null,
  h('span', { class: 'meta' }, meta));
  const li = h('li', {}, row);
  if (node.children.length) li.append(h('ul', { class: 'tree' }, node.children.map((c) => renderNode(c, activeId, onOpen))));
  return li;
}

export function historyTree(doc, docs, onOpen) {
  if (!doc?.cloning) return h('div');
  const tree = buildTree(doc, docs);
  return h('div', {},
    h('h3', {}, 'How this was made'),
    h('ul', { class: 'tree root' }, renderNode(tree, doc.id, onOpen)),
    h('p', { class: 'sub' }, 'Open the source files to follow the history further back.'));
}
