// A rotatable Cα-trace view of a protein structure, drawn on a canvas.
// Drag to rotate, scroll to zoom, double-click to reset.

import { normalise } from '../pdb.js';

const CHAIN_COLORS = ['#4f8fd6', '#e0884c', '#5aa58f', '#b07aa1', '#d6a14a', '#7a5fb0', '#c46f6f', '#69a6c9'];
const SS_COLORS = { helix: '#d6504f', sheet: '#e0b84a', coil: '#8c98a8' };

function rainbow(t) {
  const h = 240 - 240 * t; // blue at the N terminus, red at the C terminus
  return `hsl(${h} 70% 50%)`;
}

// Rotation as a 3×3 matrix, built from two angles.
const rotation = (rx, ry) => {
  const cx = Math.cos(rx); const sx = Math.sin(rx); const cy = Math.cos(ry); const sy = Math.sin(ry);
  return [[cy, 0, sy], [sx * sy, cx, -sx * cy], [-cx * sy, sx, cx * cy]];
};
const apply = (m, [x, y, z]) => [m[0][0] * x + m[0][1] * y + m[0][2] * z, m[1][0] * x + m[1][1] * y + m[1][2] * z, m[2][0] * x + m[2][1] * y + m[2][2] * z];

/**
 * @param {ReturnType<import('../pdb.js').parsePdb>} structure
 * @param {{ colour?: 'chain'|'rainbow'|'structure', highlight?: [number, number]|null, size?: number }} [options]
 */
export function proteinView(structure, { colour = 'structure', highlight = null, size = 300 } = {}) {
  const canvas = document.createElement('canvas');
  const dpr = globalThis.devicePixelRatio || 1;
  canvas.width = size * dpr;
  canvas.height = size * dpr;
  canvas.style.cssText = `width:100%;max-width:${size}px;aspect-ratio:1;display:block;margin:0 auto;border:1px solid var(--line);border-radius:8px;background:var(--paper);touch-action:none;cursor:grab`;
  canvas.setAttribute('role', 'img');
  canvas.setAttribute('aria-label', `Backbone of ${structure.name}`);
  const { centre, radius } = normalise(structure);
  let rx = -0.4; let ry = 0.6; let zoom = 1;

  const colourOf = (chainIndex, atom, i, n) => {
    if (highlight && atom.number >= highlight[0] && atom.number <= highlight[1]) return '#f2b01e';
    if (colour === 'chain') return CHAIN_COLORS[chainIndex % CHAIN_COLORS.length];
    if (colour === 'rainbow') return rainbow(i / Math.max(1, n - 1));
    return SS_COLORS[atom.ss] ?? SS_COLORS.coil;
  };

  function draw() {
    const ctx = canvas.getContext('2d');
    if (!ctx) return;
    const W = canvas.width;
    ctx.setTransform(1, 0, 0, 1, 0, 0);
    ctx.clearRect(0, 0, W, W);
    const m = rotation(rx, ry);
    const scale = (W / 2.3 / radius) * zoom;
    const segs = [];
    structure.chains.forEach((chain, ci) => {
      const pts = chain.atoms.map((a) => apply(m, [a.x - centre[0], a.y - centre[1], a.z - centre[2]]));
      for (let i = 1; i < pts.length; i++) {
        // a gap of more than 4.5 Å between neighbouring Cα atoms is a chain break
        const d = Math.hypot(chain.atoms[i].x - chain.atoms[i - 1].x, chain.atoms[i].y - chain.atoms[i - 1].y, chain.atoms[i].z - chain.atoms[i - 1].z);
        if (d > 4.5) continue;
        segs.push({ a: pts[i - 1], b: pts[i], z: (pts[i - 1][2] + pts[i][2]) / 2, color: colourOf(ci, chain.atoms[i], i, pts.length) });
      }
    });
    segs.sort((p, q) => p.z - q.z); // far to near, so nearer segments overdraw
    ctx.lineCap = 'round';
    for (const sgm of segs) {
      const near = 0.55 + 0.45 * ((sgm.z / radius + 1) / 2);
      ctx.globalAlpha = Math.max(0.25, Math.min(1, near));
      ctx.strokeStyle = sgm.color;
      ctx.lineWidth = (2.2 + 1.6 * near) * (W / 300);
      ctx.beginPath();
      ctx.moveTo(W / 2 + sgm.a[0] * scale, W / 2 - sgm.a[1] * scale);
      ctx.lineTo(W / 2 + sgm.b[0] * scale, W / 2 - sgm.b[1] * scale);
      ctx.stroke();
    }
    ctx.globalAlpha = 1;
  }

  let drag = null;
  canvas.addEventListener('pointerdown', (e) => { drag = { x: e.clientX, y: e.clientY }; canvas.style.cursor = 'grabbing'; try { canvas.setPointerCapture(e.pointerId); } catch { /* uncaptured is fine */ } });
  canvas.addEventListener('pointermove', (e) => {
    if (!drag) return;
    ry += (e.clientX - drag.x) * 0.01;
    rx += (e.clientY - drag.y) * 0.01;
    drag = { x: e.clientX, y: e.clientY };
    draw();
  });
  const end = () => { drag = null; canvas.style.cursor = 'grab'; };
  canvas.addEventListener('pointerup', end);
  canvas.addEventListener('pointercancel', end);
  canvas.addEventListener('wheel', (e) => { e.preventDefault(); zoom = Math.max(0.4, Math.min(4, zoom * (e.deltaY < 0 ? 1.1 : 0.9))); draw(); }, { passive: false });
  canvas.addEventListener('dblclick', () => { rx = -0.4; ry = 0.6; zoom = 1; draw(); });
  draw();
  return canvas;
}
