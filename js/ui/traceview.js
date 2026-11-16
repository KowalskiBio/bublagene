// Chromatogram view: the four dye traces drawn on a canvas, with the called
// bases above them. Only the visible stretch is drawn, so a 13,000-sample
// trace scrolls smoothly.

import { h, clear, fmt } from './dom.js';

// The usual convention for the four dyes. Green and red read well on either
// background; G takes the page's own ink so it stays visible in both themes.
const DYE = { A: '#22a06b', C: '#2d7ff9', G: null, T: '#e5484d' };

const BASE_ROW = 22;    // room above the curves for the called bases
const RULER_ROW = 18;   // room below for the numbers

export function traceView(traces, opts = {}) {
  const { selected = 0, onSelect, pxPerBase = 14, onZoom, highlightBase = null } = opts;
  const wrap = h('div', { class: 'tracewrap' });
  const trace = traces[selected];

  {
    const bar = h('div', { class: 'tracebar' });
    for (const [i, t] of traces.entries()) {
      bar.append(h('button', {
        class: 'btn', 'aria-pressed': String(i === selected),
        onclick: () => onSelect?.(i),
      }, t.name || `trace ${t.id}`));
    }
    bar.append(h('span', { class: 'spacer' }));
    bar.append(h('span', { class: 'sub' }, trace?.error ? '' : `${fmt(trace?.bases.length ?? 0)} bases · ${fmt(trace?.sampleCount ?? 0)} samples`));
    bar.append(h('button', { class: 'btn icon', title: 'Narrower', onclick: () => onZoom?.(Math.max(4, pxPerBase - 4)) }, '−'));
    bar.append(h('button', { class: 'btn icon', title: 'Wider', onclick: () => onZoom?.(Math.min(60, pxPerBase + 4)) }, '+'));
    wrap.append(bar);
  }

  if (!trace || trace.error) {
    wrap.append(h('p', { class: 'empty', style: { padding: '16px' } }, trace?.error ?? 'No trace to show.'));
    return wrap;
  }

  const scroller = h('div', { class: 'tracescroll' });
  const canvas = h('canvas', { class: 'tracecanvas' });
  const spacer = h('div', { class: 'tracespacer' });
  scroller.append(spacer, canvas);
  wrap.append(scroller);

  // samples per base varies; lay the trace out evenly in sample space
  const perSample = (pxPerBase * trace.bases.length) / Math.max(1, trace.sampleCount);
  const fullWidth = Math.max(1, Math.round(trace.sampleCount * perSample));
  spacer.style.width = `${fullWidth}px`;

  function draw() {
    const dpr = Math.min(2, globalThis.devicePixelRatio || 1);
    const w = scroller.clientWidth;
    const h0 = Math.max(180, scroller.clientHeight);
    canvas.width = w * dpr;
    canvas.height = h0 * dpr;
    canvas.style.width = `${w}px`;
    canvas.style.height = `${h0}px`;
    canvas.style.left = `${scroller.scrollLeft}px`;
    const ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, w, h0);

    const css = getComputedStyle(document.documentElement);
    const ink = css.getPropertyValue('--ink').trim() || '#000';
    const muted = css.getPropertyValue('--muted').trim() || '#666';
    const line = css.getPropertyValue('--line').trim() || '#ddd';
    const accent = css.getPropertyValue('--accent').trim() || '#0d6e66';

    const left = scroller.scrollLeft;
    const firstSample = Math.max(0, Math.floor(left / perSample) - 2);
    const lastSample = Math.min(trace.sampleCount, Math.ceil((left + w) / perSample) + 2);
    const x = (s) => s * perSample - left;

    const top = BASE_ROW;
    const plotH = h0 - BASE_ROW - RULER_ROW;
    let peak = 1;
    for (const c of 'ACGT') {
      const arr = trace.samples[c];
      for (let s = firstSample; s < lastSample; s++) if (arr[s] > peak) peak = arr[s];
    }
    const y = (v) => top + plotH - (v / peak) * plotH;

    // baseline
    ctx.strokeStyle = line;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(0, top + plotH + 0.5);
    ctx.lineTo(w, top + plotH + 0.5);
    ctx.stroke();

    // the four dye curves
    ctx.lineWidth = 1.25;
    for (const c of 'ACGT') {
      const arr = trace.samples[c];
      ctx.strokeStyle = DYE[c] ?? ink;
      ctx.beginPath();
      let started = false;
      for (let s = firstSample; s < lastSample; s++) {
        const px = x(s);
        const py = y(arr[s]);
        if (!started) { ctx.moveTo(px, py); started = true; } else ctx.lineTo(px, py);
      }
      ctx.stroke();
    }

    // called bases, and a tick for each
    ctx.font = '11px "IBM Plex Mono", ui-monospace, monospace';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'middle';
    for (let k = 0; k < trace.positions.length; k++) {
      const p = trace.positions[k];
      if (p < firstSample - 20 || p > lastSample + 20) continue;
      const px = x(p);
      if (px < -20 || px > w + 20) continue;
      const base = trace.bases[k] ?? 'N';
      if (highlightBase === k) {
        ctx.fillStyle = accent;
        ctx.fillRect(px - pxPerBase / 2, 2, pxPerBase, BASE_ROW - 4);
        ctx.fillStyle = css.getPropertyValue('--surface').trim() || '#fff';
      } else {
        ctx.fillStyle = DYE[base] ?? ink;
      }
      ctx.fillText(base, px, BASE_ROW / 2);
      // a faint tick down to the curves
      ctx.strokeStyle = line;
      ctx.beginPath();
      ctx.moveTo(px, BASE_ROW - 3);
      ctx.lineTo(px, BASE_ROW);
      ctx.stroke();
    }

    // numbering every ten bases
    ctx.fillStyle = muted;
    ctx.font = '10px "IBM Plex Sans", system-ui, sans-serif';
    for (let k = 9; k < trace.positions.length; k += 10) {
      const px = x(trace.positions[k]);
      if (px < -30 || px > w + 30) continue;
      ctx.fillText(String(k + 1), px, h0 - RULER_ROW / 2);
    }
  }

  let frame = null;
  const schedule = () => {
    if (frame) return;
    frame = requestAnimationFrame(() => { frame = null; draw(); });
  };
  scroller.addEventListener('scroll', schedule);
  wrap.refresh = schedule;

  // tell the caller which base was clicked
  scroller.addEventListener('click', (ev) => {
    const rect = scroller.getBoundingClientRect();
    const sample = (ev.clientX - rect.left + scroller.scrollLeft) / perSample;
    let nearest = -1;
    let bestGap = Infinity;
    for (let k = 0; k < trace.positions.length; k++) {
      const gap = Math.abs(trace.positions[k] - sample);
      if (gap < bestGap) { bestGap = gap; nearest = k; }
    }
    if (nearest >= 0 && bestGap * perSample < 30) opts.onPickBase?.(nearest, trace);
  });

  queueMicrotask(schedule);
  return wrap;
}
