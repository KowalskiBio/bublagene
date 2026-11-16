// Sequence view: both strands, feature bars, translations and cut sites.
// Only the lines in view are built, so a 35 kb sequence stays responsive.

import { h, clear, fmt } from './dom.js';
import { complement, translate, mod, rangeLength } from '../seq.js';
import { codingSequence, frameOffset, isCodingSegment, featureSpan } from '../model.js';

const NUM_W = 62;      // gutter for position numbers
const CHAR_W = 8.42;   // measured below; IBM Plex Mono at 13px
const GROUP = 10;

export function sequenceView(doc, opts = {}) {
  const { selection, onSelect, showTranslation = true, showComplement = true, sites = [], highlight } = opts;
  const L = doc.seq.length;
  // The view owns the selection while the pointer is down, so that dragging
  // never has to rebuild this element (which would drop the pointer capture).
  let sel = selection;
  const wrap = h('div', { class: 'seqwrap' });
  const inner = h('div', { style: { position: 'relative' } });
  wrap.append(inner);

  // measure one character of the real font
  const probe = h('span', { class: 'mono', style: { position: 'absolute', visibility: 'hidden', fontSize: '13px' } }, 'M'.repeat(100));
  wrap.append(probe);
  let charW = CHAR_W;

  const state = { perLine: 60, lineH: 0, first: -1, last: -1 };

  function featuresOnLine(from, to) {
    const out = [];
    for (const f of doc.features) {
      for (const seg of f.segments) {
        const len = rangeLength(seg, L);
        if (!len) continue;
        const pieces = seg.end > seg.start ? [[seg.start, seg.end]] : [[seg.start, L], [0, seg.end]];
        for (const [a, b] of pieces) {
          const s0 = Math.max(a, from);
          const s1 = Math.min(b, to);
          if (s1 > s0) out.push({ f, seg, s0, s1, startsHere: a >= from, gap: !isCodingSegment(seg) });
        }
      }
    }
    return out;
  }

  // amino-acid letters for a CDS, keyed by the base index of each codon's middle
  const aaCache = new Map();
  function aaFor(f) {
    if (aaCache.has(f.id)) return aaCache.get(f.id);
    const map = new Map();
    const coding = f.segments.filter(isCodingSegment);
    // walk the feature's own 5'->3' order collecting base positions
    const order = [];
    for (const seg of coding) {
      const len = rangeLength(seg, L);
      for (let k = 0; k < len; k++) order.push((seg.start + k) % L);
    }
    if (f.strand === -1) order.reverse();
    const prot = translate(codingSequence(doc, f));
    const off = frameOffset(f);
    for (let i = 0; i < prot.length; i++) {
      const mid = order[off + i * 3 + 1];
      if (mid !== undefined) map.set(mid, prot[i]);
    }
    aaCache.set(f.id, map);
    return map;
  }

  function buildLine(i) {
    const from = i * state.perLine;
    const to = Math.min(L, from + state.perLine);
    const top = doc.seq.slice(from, to);
    const line = h('div', { class: 'seqline', style: { position: 'absolute', left: 0, right: 0, top: `${i * state.lineH}px`, height: `${state.lineH}px` } });
    let y = 0;

    // enzyme cuts above the line
    const lineSites = sites.filter((st) => st.top > from && st.top <= to);
    if (sites.length) {
      for (const st of lineSites) {
        const cx = NUM_W + (st.top - from) * charW;
        line.append(h('div', { class: 'cutlabel', style: { left: `${cx + 2}px`, top: '0px' } }, st.enzyme.name));
        line.append(h('div', { class: 'cutmark', style: { left: `${cx}px`, top: '12px', height: '30px' } }));
      }
      y += 14;
    }

    line.append(h('div', { class: 'num', style: { top: `${y + 2}px` } }, fmt(from + 1)));

    // selection highlight
    if (sel && rangeLength(sel, L) > 0) {
      const pieces = sel.end > sel.start ? [[sel.start, sel.end]] : [[sel.start, L], [0, sel.end]];
      for (const [a, b] of pieces) {
        const s0 = Math.max(a, from);
        const s1 = Math.min(b, to);
        if (s1 > s0) {
          line.append(h('div', {
            class: 'sel',
            style: { left: `${NUM_W + (s0 - from) * charW}px`, width: `${(s1 - s0) * charW}px`, top: `${y}px`, height: `${showComplement ? 34 : 17}px` },
          }));
        }
      }
    }

    line.append(h('div', { class: 'bases', style: { left: `${NUM_W}px`, top: `${y}px` } }, top));
    y += 17;
    if (showComplement) {
      line.append(h('div', { class: 'bases bottom', style: { left: `${NUM_W}px`, top: `${y}px` } }, complement(top)));
      y += 17;
    }

    // translations
    if (showTranslation) {
      for (const f of doc.features.filter((x) => x.type === 'CDS')) {
        const map = aaFor(f);
        let text = '';
        let any = false;
        for (let p = from; p < to; p++) {
          const aa = map.get(p);
          if (aa) { text += aa; any = true; } else text += ' ';
        }
        if (any) {
          line.append(h('div', { class: 'aa', style: { left: `${NUM_W}px`, top: `${y}px` }, title: f.name }, text));
          y += 14;
        }
      }
    }

    // feature bars
    for (const { f, s0, s1, startsHere, gap } of featuresOnLine(from, to)) {
      line.append(h('div', {
        class: 'fbar',
        title: f.name,
        dataset: { feature: f.id },
        style: {
          left: `${NUM_W + (s0 - from) * charW}px`,
          width: `${Math.max(2, (s1 - s0) * charW)}px`,
          top: `${y}px`,
          background: gap ? 'transparent' : f.color,
          borderTop: gap ? '2px dotted var(--muted)' : 'none',
          opacity: highlight && highlight !== f.id ? 0.4 : 1,
        },
      }, gap ? '' : (startsHere ? `${f.strand === -1 ? '◀ ' : ''}${f.name}${f.strand === 1 ? ' ▶' : ''}` : '')));
      y += 15;
    }
    return line;
  }

  // Every line is the same height, so measure the busiest one rather than
  // reserving room for every feature in the document.
  function lineHeight() {
    const total = Math.max(1, Math.ceil(L / state.perLine));
    const bars = new Int32Array(total);
    const aas = new Int32Array(total);
    const mark = (a, b, arr) => {
      const first = Math.floor(a / state.perLine);
      const last = Math.floor(Math.max(a, b - 1) / state.perLine);
      for (let i = first; i <= last && i < total; i++) arr[i]++;
    };
    const pieces = (r) => (r.end > r.start ? [[r.start, r.end]] : [[r.start, L], [0, r.end]]);
    for (const f of doc.features) {
      for (const seg of f.segments) {
        if (!rangeLength(seg, L)) continue;
        for (const [a, b] of pieces(seg)) mark(a, b, bars);
      }
      if (showTranslation && f.type === 'CDS') for (const [a, b] of pieces(featureSpan(f))) mark(a, b, aas);
    }
    let maxBars = 1;
    let maxAA = 0;
    for (let i = 0; i < total; i++) {
      if (bars[i] > maxBars) maxBars = bars[i];
      if (aas[i] > maxAA) maxAA = aas[i];
    }
    const strands = 1 + (showComplement ? 1 : 0);
    return Math.round((sites.length ? 14 : 0) + strands * 17 + maxAA * 14 + maxBars * 15 + 10);
  }

  function layout() {
    charW = probe.getBoundingClientRect().width / 100 || CHAR_W;
    const avail = wrap.clientWidth - NUM_W - 24;
    const perLine = Math.max(GROUP, Math.floor(avail / charW / GROUP) * GROUP);
    const lineH = lineHeight();
    if (perLine === state.perLine && lineH === state.lineH && state.first >= 0) return false;
    state.perLine = perLine;
    state.lineH = lineH;
    state.first = -1;
    inner.style.height = `${Math.ceil(L / perLine) * lineH}px`;
    return true;
  }

  function renderWindow() {
    const scroller = wrap.parentElement;
    if (!scroller) return;
    const top = scroller.scrollTop - wrap.offsetTop;
    const height = scroller.clientHeight;
    const total = Math.ceil(L / state.perLine);
    const first = Math.max(0, Math.floor(top / state.lineH) - 3);
    const last = Math.min(total - 1, Math.ceil((top + height) / state.lineH) + 3);
    if (first === state.first && last === state.last) return;
    state.first = first;
    state.last = last;
    clear(inner);
    for (let i = first; i <= last; i++) inner.append(buildLine(i));
  }

  // click and drag to select
  let anchor = null;
  const posAt = (ev) => {
    const r = inner.getBoundingClientRect();
    const line = Math.floor((ev.clientY - r.top) / state.lineH);
    const col = Math.round((ev.clientX - r.left - NUM_W) / charW);
    if (line < 0) return 0;
    const clamped = Math.max(0, Math.min(state.perLine, col));
    return Math.max(0, Math.min(L, line * state.perLine + clamped));
  };
  const rangeFrom = (p) => (p === anchor ? null : { start: Math.min(anchor, p), end: Math.max(anchor, p) });

  wrap.addEventListener('pointerdown', (ev) => {
    if (ev.button !== 0) return;
    anchor = posAt(ev);
    wrap.setPointerCapture(ev.pointerId);
    ev.preventDefault();
  });
  wrap.addEventListener('pointermove', (ev) => {
    if (anchor === null) return;
    const next = rangeFrom(posAt(ev));
    if (next?.start === sel?.start && next?.end === sel?.end) return;
    wrap.setSelection(next);
    onSelect?.(next, { live: true });
  });
  const finish = (ev) => {
    if (anchor === null) return;
    const next = rangeFrom(posAt(ev));
    anchor = null;
    wrap.setSelection(next);
    onSelect?.(next, { live: false });
  };
  wrap.addEventListener('pointerup', finish);
  wrap.addEventListener('pointercancel', finish);

  wrap.posFromEvent = posAt;
  // repaint the visible lines in place, keeping this element alive
  wrap.setSelection = (next) => {
    sel = next;
    state.first = -1;
    renderWindow();
  };
  wrap.refresh = () => { if (layout()) { state.first = -1; } renderWindow(); };
  wrap.onScroll = renderWindow;
  queueMicrotask(() => wrap.refresh());
  return wrap;
}
