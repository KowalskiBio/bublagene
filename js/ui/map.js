// Circular and linear sequence maps, drawn as SVG.

import { s, niceStep, fmt } from './dom.js';
import { featureSpan, featureLength } from '../model.js';
import { rangeLength, mod } from '../seq.js';

const TAU = Math.PI * 2;

// base index -> angle in radians, origin at 12 o'clock, running clockwise
const angleOf = (i, L) => (i / L) * TAU - Math.PI / 2;
const pt = (cx, cy, r, a) => [cx + r * Math.cos(a), cy + r * Math.sin(a)];

// An annular wedge from a0 to a1, optionally pointed at one end.
function arcPath(cx, cy, a0, a1, rIn, rOut, arrow = 0) {
  const sweep = a1 - a0;
  const big = sweep > Math.PI ? 1 : 0;
  const rMid = (rIn + rOut) / 2;
  const tipA = arrow > 0 ? a1 : a0;
  const bodyEnd = arrow > 0 ? a1 - arrow : a0 + arrow;
  if (!arrow) {
    const [x0, y0] = pt(cx, cy, rOut, a0);
    const [x1, y1] = pt(cx, cy, rOut, a1);
    const [x2, y2] = pt(cx, cy, rIn, a1);
    const [x3, y3] = pt(cx, cy, rIn, a0);
    return `M${x0} ${y0}A${rOut} ${rOut} 0 ${big} 1 ${x1} ${y1}L${x2} ${y2}A${rIn} ${rIn} 0 ${big} 0 ${x3} ${y3}Z`;
  }
  const from = arrow > 0 ? a0 : bodyEnd;
  const to = arrow > 0 ? bodyEnd : a1;
  const b2 = Math.abs(to - from) > Math.PI ? 1 : 0;
  const [x0, y0] = pt(cx, cy, rOut, from);
  const [x1, y1] = pt(cx, cy, rOut, to);
  const [tx, ty] = pt(cx, cy, rMid, tipA);
  const [x2, y2] = pt(cx, cy, rIn, to);
  const [x3, y3] = pt(cx, cy, rIn, from);
  if (arrow > 0) {
    return `M${x0} ${y0}A${rOut} ${rOut} 0 ${b2} 1 ${x1} ${y1}L${tx} ${ty}L${x2} ${y2}A${rIn} ${rIn} 0 ${b2} 0 ${x3} ${y3}Z`;
  }
  return `M${tx} ${ty}L${x0} ${y0}A${rOut} ${rOut} 0 ${b2} 1 ${x1} ${y1}L${x2} ${y2}A${rIn} ${rIn} 0 ${b2} 0 ${x3} ${y3}Z`;
}

// Keep labels from overlapping: push them apart along y, per side.
function spread(items, minGap, lo, hi) {
  items.sort((a, b) => a.y - b.y);
  for (let i = 1; i < items.length; i++) {
    if (items[i].y - items[i - 1].y < minGap) items[i].y = items[i - 1].y + minGap;
  }
  const over = items.length ? items[items.length - 1].y - hi : 0;
  if (over > 0) for (const it of items) it.y -= over;
  for (let i = items.length - 2; i >= 0; i--) {
    if (items[i + 1].y - items[i].y < minGap) items[i].y = items[i + 1].y - minGap;
  }
  if (items.length && items[0].y < lo) {
    const d = lo - items[0].y;
    for (const it of items) it.y += d;
  }
  return items;
}

// Stack features that overlap into concentric lanes, so none is hidden.
function assignLanes(features, L) {
  const taken = [];
  const out = new Map();
  const ordered = [...features].sort((a, b) => featureSpan(a).start - featureSpan(b).start);
  for (const f of ordered) {
    const span = featureSpan(f);
    const a0 = angleOf(span.start, L);
    const a1 = a0 + (rangeLength(span, L) / L) * TAU;
    let lane = 0;
    while (lane < 6 && (taken[lane] ??= []).some(([b0, b1]) => a0 < b1 - 1e-4 && b0 < a1 - 1e-4)) lane++;
    (taken[lane] ??= []).push([a0, a1]);
    out.set(f.id, Math.min(lane, 5));
  }
  return { lanes: out, depth: Math.max(1, taken.length) };
}

export function circularMap(doc, opts = {}) {
  const { selection, sites = [], highlight, onPick, onPickSite, overlays = [] } = opts;
  const features = overlays.length ? [...doc.features, ...overlays] : doc.features;
  const L = doc.seq.length || 1;
  const W = 680;
  const cx = W / 2;
  const cy = W / 2;
  const BAND = 15;          // thickness of a feature arrow
  const LANE = BAND + 4;
  const { lanes, depth } = assignLanes(features, L);
  const R = 150;            // the backbone, with features stacked outward from it
  const outerR = R + depth * LANE;
  const root = s('svg', { viewBox: `0 0 ${W} ${W}`, width: W, height: W, role: 'img', 'aria-label': `Map of ${doc.name}` });

  // backbone
  root.append(s('circle', { cx, cy, r: R, fill: 'none', stroke: 'var(--line)', 'stroke-width': 3 }));
  if (!doc.circular) {
    root.append(s('path', { d: `M${cx} ${cy - R - 6}L${cx} ${cy - R + 6}`, stroke: 'var(--paper)', 'stroke-width': 7 }));
  }

  // selection arc, on the backbone itself
  if (selection && rangeLength(selection, L) > 0) {
    const a0 = angleOf(selection.start, L);
    const a1 = a0 + (rangeLength(selection, L) / L) * TAU;
    root.append(s('path', { d: arcPath(cx, cy, a0, a1, R - 7, R + 7), fill: 'var(--select)', opacity: 0.6 }));
  }

  // ruler, inside the backbone so it never runs into a feature
  const step = niceStep(L, 8);
  for (let p = 0; p < L; p += step) {
    const a = angleOf(p, L);
    const [x0, y0] = pt(cx, cy, R - 3, a);
    const [x1, y1] = pt(cx, cy, R - 11, a);
    root.append(s('line', { x1: x0, y1: y0, x2: x1, y2: y1, stroke: 'var(--line)', 'stroke-width': 1 }));
    const [tx, ty] = pt(cx, cy, R - 23, a);
    root.append(s('text', {
      x: tx, y: ty, 'font-size': 9.5, fill: 'var(--muted)', 'text-anchor': 'middle', 'dominant-baseline': 'middle',
    }, p === 0 ? '1' : fmt(p)));
  }

  // features, as arrows on their own lane
  const labels = [];
  for (const f of features) {
    const lane = lanes.get(f.id) ?? 0;
    const rMid = R + (lane + 0.5) * LANE;
    const rIn = rMid - BAND / 2;
    const rOut = rMid + BAND / 2;
    const span = featureSpan(f);
    for (const seg of f.segments) {
      const segLen = rangeLength(seg, L);
      if (!segLen) continue;
      const a0 = angleOf(seg.start, L);
      const a1 = a0 + (segLen / L) * TAU;
      const isTip = f.strand === 1 ? seg === f.segments[f.segments.length - 1] : seg === f.segments[0];
      const arrowSpan = Math.min((a1 - a0) * 0.4, 14 / rMid);
      const arrow = f.strand && isTip ? (f.strand === 1 ? arrowSpan : -arrowSpan) : 0;
      const gap = seg.extra?.type === 'gap';
      const path = s('path', {
        d: gap ? arcPath(cx, cy, a0, a1, rMid - 1, rMid + 1) : arcPath(cx, cy, a0, a1, rIn, rOut, arrow),
        fill: f.orf ? 'none' : gap ? 'var(--muted)' : (seg.color && seg.color !== 'noColor' ? seg.color : f.color),
        opacity: highlight && highlight !== f.id ? 0.3 : 1,
        stroke: highlight === f.id ? 'var(--ink)' : f.orf ? 'var(--muted)' : 'none',
        'stroke-dasharray': f.orf ? '3 2' : null,
        'stroke-width': f.orf ? 1.2 : 1.5,
        style: 'cursor:pointer',
        'data-feature': f.id,
      });
      path.append(s('title', {}, `${f.name} · ${fmt(featureLength(f, L))} bp`));
      if (onPick) path.addEventListener('click', (e) => { e.stopPropagation(); onPick(f); });
      root.append(path);
    }
    const mid = angleOf(span.start, L) + (rangeLength(span, L) / L) * TAU / 2;
    if (!f.orf) labels.push({ a: mid, text: f.name, id: f.id, from: R + (lane + 1) * LANE });
  }

  // enzyme ticks, clear of every feature lane
  for (const site of sites) {
    const a = angleOf(site.top, L);
    const [x0, y0] = pt(cx, cy, outerR + 3, a);
    const [x1, y1] = pt(cx, cy, outerR + 11, a);
    root.append(s('line', { x1: x0, y1: y0, x2: x1, y2: y1, stroke: 'var(--muted)', 'stroke-width': 1 }));
    labels.push({ a, text: `${site.enzyme.name} (${fmt(site.top)})`, site, kind: 'site', from: outerR + 11 });
  }

  // Every leader line starts on one circle, so they fan out without crossing.
  const anchorR = outerR + 14;
  const left = [];
  const right = [];
  for (const lab of labels) {
    const [ax, ay] = pt(cx, cy, anchorR, lab.a);
    const [fx, fy] = pt(cx, cy, lab.from, lab.a);
    (Math.cos(lab.a) < 0 ? left : right).push({ ...lab, x: ax, y: ay, anchorX: ax, anchorY: ay, footX: fx, footY: fy });
  }
  spread(left, 15, 14, W - 14);
  spread(right, 15, 14, W - 14);
  for (const side of [left, right]) {
    const isLeft = side === left;
    for (const lab of side) {
      const x = isLeft ? 14 : W - 14;
      const elbow = isLeft ? x + 30 : x - 30;
      root.append(s('polyline', {
        points: `${lab.footX},${lab.footY} ${lab.anchorX},${lab.anchorY} ${elbow},${lab.y} ${isLeft ? x + 3 : x - 3},${lab.y}`,
        fill: 'none', stroke: 'var(--line)', 'stroke-width': 1,
      }));
      const t = s('text', {
        x, y: lab.y, 'font-size': 11.5, fill: lab.kind === 'site' ? 'var(--muted)' : 'var(--ink)',
        'text-anchor': isLeft ? 'start' : 'end', 'dominant-baseline': 'middle',
        style: 'cursor:pointer',
        'font-weight': lab.kind === 'site' ? 400 : 600,
      }, lab.text);
      if (lab.kind === 'site' && onPickSite) t.addEventListener('click', () => onPickSite(lab.site));
      else if (onPick) t.addEventListener('click', () => onPick(features.find((f) => f.id === lab.id)));
      root.append(t);
    }
  }

  // centre caption
  root.append(s('text', { x: cx, y: cy - 7, 'text-anchor': 'middle', 'font-size': 16, 'font-weight': 600, fill: 'var(--ink)' }, doc.name));
  root.append(s('text', { x: cx, y: cy + 14, 'text-anchor': 'middle', 'font-size': 11.5, fill: 'var(--muted)' }, `${fmt(L)} bp · ${doc.circular ? 'circular' : 'linear'}`));
  return root;
}

export function linearMap(doc, opts = {}) {
  const { selection, sites = [], highlight, onPick, onPickSite, overlays = [] } = opts;
  const features = overlays.length ? [...doc.features, ...overlays] : doc.features;
  const L = doc.seq.length || 1;
  const W = 900;
  const padX = 52;
  const trackY = 150;
  const innerW = W - padX * 2;
  const x = (p) => padX + (p / L) * innerW;
  const rows = Math.max(1, Math.min(4, features.length));
  const H = trackY + 70 + rows * 6;
  const root = s('svg', { viewBox: `0 0 ${W} ${H}`, width: W, height: H, role: 'img', 'aria-label': `Map of ${doc.name}` });

  if (selection && rangeLength(selection, L) > 0) {
    const segs = selection.end > selection.start
      ? [[selection.start, selection.end]]
      : [[selection.start, L], [0, selection.end]];
    for (const [a, b] of segs) {
      root.append(s('rect', { x: x(a), y: trackY - 14, width: Math.max(1, x(b) - x(a)), height: 28, fill: 'var(--select)', opacity: 0.55 }));
    }
  }

  root.append(s('line', { x1: padX, y1: trackY, x2: W - padX, y2: trackY, stroke: 'var(--line)', 'stroke-width': 7, 'stroke-linecap': doc.circular ? 'butt' : 'round' }));

  // features on stacked lanes so overlapping ones stay readable
  const lanes = [];
  const placed = [];
  for (const f of features) {
    const span = featureSpan(f);
    const x0 = x(span.start);
    const x1 = x(span.end > span.start ? span.end : L);
    let lane = 0;
    while (lanes[lane] !== undefined && lanes[lane] > x0 - 4) lane++;
    lanes[lane] = Math.max(x1, x0 + 40);
    placed.push({ f, lane, x0, x1 });
  }
  for (const { f, lane, x0, x1 } of placed) {
    const y = trackY + (f.strand === -1 ? 20 : -34) + (f.strand === -1 ? lane * 20 : -lane * 20);
    const w = Math.max(3, x1 - x0);
    const g = s('g', { style: 'cursor:pointer', opacity: highlight && highlight !== f.id ? 0.4 : 1, 'data-feature': f.id });
    const tip = Math.min(8, w / 2);
    const d = f.strand === -1
      ? `M${x0 + w} ${y}L${x0 + tip} ${y}L${x0} ${y + 7}L${x0 + tip} ${y + 14}L${x0 + w} ${y + 14}Z`
      : `M${x0} ${y}L${x0 + w - tip} ${y}L${x0 + w} ${y + 7}L${x0 + w - tip} ${y + 14}L${x0} ${y + 14}Z`;
    g.append(s('path', { d, fill: f.orf ? 'none' : f.color, stroke: highlight === f.id ? 'var(--ink)' : f.orf ? 'var(--muted)' : 'none', 'stroke-dasharray': f.orf ? '3 2' : null, 'stroke-width': f.orf ? 1.2 : 1.5 }));
    g.append(s('text', { x: x0 + w / 2, y: y + (f.strand === -1 ? 26 : -4), 'font-size': 11, 'text-anchor': 'middle', fill: 'var(--ink)' }, f.name));
    g.append(s('title', {}, `${f.name} · ${fmt(featureLength(f, L))} bp`));
    if (onPick) g.addEventListener('click', () => onPick(f));
    root.append(g);
  }

  // ruler
  const step = niceStep(L, 9);
  for (let p = 0; p <= L; p += step) {
    root.append(s('line', { x1: x(p), y1: trackY + 8, x2: x(p), y2: trackY + 14, stroke: 'var(--muted)' }));
    root.append(s('text', { x: x(p), y: trackY + 27, 'font-size': 10, fill: 'var(--muted)', 'text-anchor': 'middle' }, p === 0 ? '1' : fmt(p)));
  }

  // enzyme sites above
  const used = [];
  for (const site of sites) {
    const sx = x(site.top);
    let row = 0;
    while (used[row] !== undefined && used[row] > sx - 46) row++;
    used[row] = sx;
    const y = 56 - row * 15;
    root.append(s('line', { x1: sx, y1: y + 6, x2: sx, y2: trackY - 16, stroke: 'var(--line)', 'stroke-dasharray': '2 2' }));
    const t = s('text', { x: sx, y, 'font-size': 10.5, fill: 'var(--muted)', 'text-anchor': 'middle', style: 'cursor:pointer' }, site.enzyme.name);
    if (onPickSite) t.addEventListener('click', () => onPickSite(site));
    root.append(t);
  }
  return root;
}
