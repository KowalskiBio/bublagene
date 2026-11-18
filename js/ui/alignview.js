// @ts-nocheck — canvas and DOM drawing.
// Alignment view: the reference as a map, the aligned sequences beneath it
// column by column, and the chromatogram of the selected one under that, all
// sharing one horizontal scroll so a peak sits under the base it called.

import { complement } from '../seq.js';
import { clear, fmt, h, s } from './dom.js';

const GUTTER = 120;      // room for row names, left of the first base
const RULER = 16;
const ROW = 18;
const TRACE_H = 170;
const BASE_ROW = 20;     // called bases above the curves
const DYE = { A: '#22a06b', C: '#2d7ff9', G: null, T: '#e5484d' };

// Pack the reference's features into rows so overlapping ones do not hide each other.
function packFeatures(features) {
  const items = features
    .flatMap((f) => f.segments.map((seg) => ({ f, start: seg.start, end: seg.end })))
    .sort((a, b) => a.start - b.start);
  const lanes = [];
  for (const it of items) {
    let lane = lanes.findIndex((end) => end <= it.start);
    if (lane < 0) { lane = lanes.length; lanes.push(0); }
    lanes[lane] = it.end;
    it.lane = lane;
  }
  return { items, depth: Math.max(1, lanes.length) };
}

function overview(doc, rows, selected, scroller, geometry, onJump) {
  const L = doc.seq.length;
  const { items, depth } = packFeatures(doc.features);
  const BAR = 11;
  const H = 8 + depth * (BAR + 2) + rows.length * 5 + 6;
  const X = 1000;
  const svg = s('svg', { viewBox: `0 0 ${X} ${H}`, preserveAspectRatio: 'none', class: 'alignmapsvg', role: 'img', 'aria-label': `Map of ${doc.name}`, style: `height:${H}px` });
  svg.append(s('rect', { x: 0, y: 0, width: X, height: H, fill: 'var(--surface)' }));
  for (const it of items) {
    const x = (it.start / L) * X;
    const w = Math.max(1.5, ((it.end - it.start) / L) * X);
    const y = 4 + it.lane * (BAR + 2);
    const r = s('rect', { x, y, width: w, height: BAR, rx: 2, fill: it.f.color, opacity: 0.9 });
    r.append(s('title', {}, `${it.f.name} · ${it.start + 1}–${it.end}`));
    svg.append(r);
  }
  rows.forEach((row, i) => {
    const x = (row.start / L) * X;
    const w = Math.max(1.5, ((row.end - row.start) / L) * X);
    svg.append(s('rect', { x, y: 6 + depth * (BAR + 2) + i * 5, width: w, height: 3, rx: 1.5, fill: i === selected ? 'var(--accent)' : 'var(--muted)', opacity: i === selected ? 1 : 0.55 }));
  });
  const win = s('rect', { y: 0, height: H, fill: 'var(--accent)', opacity: 0.16, stroke: 'var(--accent)', 'stroke-width': 1, 'vector-effect': 'non-scaling-stroke' });
  svg.append(win);
  const place = () => {
    const per = geometry.cw;
    const first = Math.max(0, scroller.scrollLeft / per);
    const last = Math.min(L, (scroller.scrollLeft + scroller.clientWidth - GUTTER) / per);
    win.setAttribute('x', String((first / L) * X));
    win.setAttribute('width', String(Math.max(2, ((last - first) / L) * X)));
  };
  svg.addEventListener('click', (e) => {
    const rect = svg.getBoundingClientRect();
    onJump(((e.clientX - rect.left) / rect.width) * L);
  });
  svg.place = place;
  return svg;
}

export function alignmentView(doc, rows, traces, opts = {}) {
  const { selected = 0, cw = 10, onSelect, onZoom, onPickColumn, scrollTo = null, onScroll } = opts;
  const L = doc.seq.length;
  const ref = doc.seq.toUpperCase();
  const row = rows[Math.min(selected, rows.length - 1)];
  const trace = row?.traceId ? traces.find((t) => t.id === row.traceId) : null;
  const wrap = h('div', { class: 'alignwrap' });

  const bar = h('div', { class: 'alignbar' });
  rows.forEach((r, i) => bar.append(h('button', {
    class: 'btn', 'aria-pressed': String(i === selected), onclick: () => onSelect?.(i),
    title: `${fmt(r.start + 1)}–${fmt(r.end)}${r.identity != null ? ` · ${r.identity.toFixed(1)}% identity` : ''}`,
  }, `${r.reverse ? '◀' : '▶'} ${r.label}`)));
  bar.append(h('span', { class: 'spacer' }));
  if (row) bar.append(h('span', { class: 'sub' }, `${fmt(row.end - row.start)} columns${row.identity != null ? ` · ${row.identity.toFixed(1)}% identity` : ''}${row.reverse ? ' · reverse strand' : ''}`));
  bar.append(h('button', { class: 'btn icon', title: 'Narrower', onclick: () => onZoom?.(Math.max(4, cw - 2)) }, '−'));
  bar.append(h('button', { class: 'btn icon', title: 'Wider', onclick: () => onZoom?.(Math.min(24, cw + 2)) }, '+'));
  wrap.append(bar);

  const geometry = { cw };
  const scroller = h('div', { class: 'alignscroll' });
  const canvas = h('canvas', { class: 'aligncanvas' });
  const spacer = h('div', { class: 'alignspacer', style: { width: `${GUTTER + L * cw + 40}px` } });
  scroller.append(spacer, canvas);

  const map = overview(doc, rows, selected, scroller, geometry, (pos) => {
    scroller.scrollLeft = Math.max(0, pos * cw - (scroller.clientWidth - GUTTER) / 2);
  });
  wrap.append(map, scroller);

  const rowsH = RULER + ROW + rows.length * ROW + 10;
  const totalH = rowsH + TRACE_H;
  scroller.style.height = `${totalH + 18}px`;

  // x position of each called base of the selected trace, in the same columns as the alignment
  let xs = null;
  let peak = 1;
  if (trace?.refMap && trace.positions?.length) {
    const n = trace.bases.length;
    xs = new Array(n).fill(null);
    for (let i = 0; i < n; i++) if (trace.refMap[i] != null) xs[i] = (trace.refMap[i] + 0.5) * cw;
    // insertions share a column with a neighbour; give them a little room of their own
    for (let i = 1; i < n; i++) {
      if (xs[i] != null && xs[i - 1] != null && Math.abs(xs[i] - xs[i - 1]) < cw * 0.2) {
        xs[i] = xs[i - 1] + (xs[i] >= xs[i - 1] ? 1 : -1) * cw * 0.3;
      }
    }
    for (const c of 'ACGT') for (const v of trace.samples[c]) if (v > peak) peak = v;
  }

  function draw() {
    const dpr = Math.min(2, globalThis.devicePixelRatio || 1);
    const w = scroller.clientWidth;
    const hh = totalH;
    canvas.width = w * dpr;
    canvas.height = hh * dpr;
    canvas.style.width = `${w}px`;
    canvas.style.height = `${hh}px`;
    canvas.style.left = `${scroller.scrollLeft}px`;
    const ctx = canvas.getContext('2d');
    if (!ctx) return;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, hh);
    const css = getComputedStyle(document.documentElement);
    const ink = css.getPropertyValue('--ink').trim() || '#000';
    const muted = css.getPropertyValue('--muted').trim() || '#666';
    const line = css.getPropertyValue('--line').trim() || '#ddd';
    const accent = css.getPropertyValue('--accent').trim() || '#0d6e66';
    const danger = css.getPropertyValue('--danger').trim() || '#b00';
    const mismatch = css.getPropertyValue('--mismatch').trim() || '#b45309';
    const sunk = css.getPropertyValue('--sunk').trim() || '#eee';
    const left = scroller.scrollLeft;
    const xOf = (col) => GUTTER + col * cw - left;
    const first = Math.max(0, Math.floor((left) / cw) - 1);
    const last = Math.min(L, Math.ceil((left + w - GUTTER) / cw) + 1);
    const colour = (b) => DYE[b] ?? ink;
    const fontPx = Math.max(9, Math.min(13, cw + 1));
    const letters = cw >= 7;

    // clip the data area so nothing is drawn over the names
    ctx.save();
    ctx.beginPath();
    ctx.rect(GUTTER, 0, w - GUTTER, hh);
    ctx.clip();

    // ruler and the reference's own bases
    ctx.fillStyle = muted;
    ctx.font = '10px "IBM Plex Sans", system-ui, sans-serif';
    ctx.textBaseline = 'middle';
    ctx.textAlign = 'center';
    const step = cw >= 12 ? 10 : cw >= 7 ? 10 : 20;
    for (let p = Math.ceil((first + 1) / step) * step; p <= last; p += step) {
      ctx.fillText(fmt(p), xOf(p - 1) + cw / 2, RULER / 2);
      ctx.fillRect(xOf(p - 1) + cw / 2 - 0.5, RULER - 3, 1, 3);
    }
    ctx.font = `${fontPx}px "IBM Plex Mono", ui-monospace, monospace`;
    for (let c = first; c < last; c++) {
      const x = xOf(c) + cw / 2;
      if (letters) { ctx.fillStyle = colour(ref[c]); ctx.fillText(ref[c], x, RULER + ROW / 2); }
      else { ctx.fillStyle = colour(ref[c]); ctx.globalAlpha = 0.7; ctx.fillRect(xOf(c), RULER + 4, Math.max(1, cw - 1), ROW - 8); ctx.globalAlpha = 1; }
    }

    // each aligned sequence, mismatches filled in
    rows.forEach((r, i) => {
      const y = RULER + ROW + i * ROW;
      if (i === selected) { ctx.fillStyle = sunk; ctx.fillRect(GUTTER, y, w - GUTTER, ROW); }
      for (let c = Math.max(first, r.start); c < Math.min(last, r.end); c++) {
        const ch = r.chars[c - r.start];
        const x = xOf(c);
        if (ch === '-') {
          ctx.fillStyle = danger;
          ctx.fillRect(x + 1, y + ROW / 2 - 1, Math.max(2, cw - 2), 2);
        } else if (ch === ref[c]) {
          ctx.globalAlpha = i === selected ? 0.95 : 0.6;
          ctx.fillStyle = colour(ch);
          if (letters) ctx.fillText(ch, x + cw / 2, y + ROW / 2); else ctx.fillRect(x, y + 4, Math.max(1, cw - 1), ROW - 8);
          ctx.globalAlpha = 1;
        } else {
          ctx.fillStyle = mismatch;
          ctx.fillRect(x, y + 1, Math.max(2, cw - 1), ROW - 2);
          ctx.fillStyle = '#fff';
          if (letters) ctx.fillText(ch, x + cw / 2, y + ROW / 2);
        }
      }
      for (const [at] of r.insertions) {
        if (at < first || at > last) continue;
        ctx.fillStyle = accent;
        ctx.fillRect(xOf(at) - 1, y + 1, 2, ROW - 2);
      }
    });

    // the chromatogram of the selected sequence, stretched onto the same columns
    const top = rowsH + BASE_ROW;
    const plotH = TRACE_H - BASE_ROW - 14;
    ctx.fillStyle = line;
    ctx.fillRect(GUTTER, rowsH, w - GUTTER, 1);
    if (xs) {
      const y = (v) => top + plotH - (v / peak) * plotH;
      const flip = !!trace.reverse;
      const shade = (c) => DYE[flip ? complement(c) : c] ?? ink;
      ctx.lineWidth = 1.25;
      for (const c of 'ACGT') {
        ctx.strokeStyle = shade(c);
        ctx.beginPath();
        const arr = trace.samples[c];
        let open = false;
        for (let i = 0; i + 1 < xs.length; i++) {
          const x0 = xs[i];
          const x1 = xs[i + 1];
          if (x0 == null || x1 == null) { open = false; continue; }
          const lo = Math.min(x0, x1) - left + GUTTER;
          const hi = Math.max(x0, x1) - left + GUTTER;
          if (hi < GUTTER - 4 || lo > w + 4) { open = false; continue; }
          const p0 = trace.positions[i];
          const p1 = trace.positions[i + 1];
          const span = Math.max(1, p1 - p0);
          const stepS = Math.max(1, Math.floor(span / 24));
          for (let sIdx = p0; sIdx <= p1; sIdx += stepS) {
            const px = GUTTER + x0 + ((sIdx - p0) / span) * (x1 - x0) - left;
            const py = y(arr[sIdx] ?? 0);
            if (!open) { ctx.moveTo(px, py); open = true; } else ctx.lineTo(px, py);
          }
        }
        ctx.stroke();
      }
      ctx.font = `${fontPx}px "IBM Plex Mono", ui-monospace, monospace`;
      ctx.textAlign = 'center';
      for (let i = 0; i < xs.length; i++) {
        if (xs[i] == null) continue;
        const px = GUTTER + xs[i] - left;
        if (px < GUTTER - 10 || px > w + 10) continue;
        const b = trace.bases[i] ?? 'N';
        const shown = flip ? complement(b) : b;
        ctx.fillStyle = DYE[shown] ?? ink;
        if (letters) ctx.fillText(shown, px, rowsH + BASE_ROW / 2);
        else ctx.fillRect(px - 1, rowsH + 5, 2, 8);
      }
    } else {
      ctx.fillStyle = muted;
      ctx.font = '12px "IBM Plex Sans", system-ui, sans-serif';
      ctx.textAlign = 'left';
      ctx.fillText(row ? 'No chromatogram for this sequence.' : '', GUTTER + 12, rowsH + TRACE_H / 2);
    }
    ctx.restore();

    // names, over a plain background so scrolled bases never run under them
    const paper = css.getPropertyValue('--surface').trim() || '#fff';
    ctx.fillStyle = paper;
    ctx.fillRect(0, 0, GUTTER - 4, hh);
    ctx.font = '11px "IBM Plex Sans", system-ui, sans-serif';
    ctx.textAlign = 'left';
    ctx.textBaseline = 'middle';
    const name = (t, max) => (t.length > max ? `${t.slice(0, max - 1)}…` : t);
    ctx.fillStyle = muted;
    ctx.fillText(name(doc.name, 17), 8, RULER + ROW / 2);
    rows.forEach((r, i) => {
      ctx.fillStyle = i === selected ? ink : muted;
      ctx.font = `${i === selected ? 600 : 400} 11px "IBM Plex Sans", system-ui, sans-serif`;
      ctx.fillText(`${r.reverse ? '◀' : '▶'} ${name(r.label, 15)}`, 8, RULER + ROW + i * ROW + ROW / 2);
    });
    ctx.fillStyle = muted;
    ctx.font = '11px "IBM Plex Sans", system-ui, sans-serif';
    ctx.fillText(trace ? name(trace.name || 'trace', 17) : '', 8, rowsH + TRACE_H / 2);
    map.place();
  }

  let frame = null;
  const schedule = () => { if (frame) return; frame = requestAnimationFrame(() => { frame = null; draw(); }); };
  scroller.addEventListener('scroll', () => { onScroll?.(scroller.scrollLeft); schedule(); });
  wrap.refresh = schedule;

  scroller.addEventListener('click', (e) => {
    const rect = scroller.getBoundingClientRect();
    const x = e.clientX - rect.left;
    if (x < GUTTER) return;
    const col = Math.floor((x - GUTTER + scroller.scrollLeft) / cw);
    if (col >= 0 && col < L) onPickColumn?.(col);
  });

  queueMicrotask(() => {
    if (scrollTo != null) scroller.scrollLeft = scrollTo;
    else if (row) scroller.scrollLeft = Math.max(0, row.start * cw - 40);
    schedule();
  });
  return wrap;
}
