// Sequence view: both strands, feature bars, translations and cut sites.
// Only the lines in view are built, so a 35 kb sequence stays responsive.

import { h, clear, fmt } from './dom.js';
import { complement, translate, revcomp, mod, rangeLength } from '../seq.js';
import { codingSequence, frameOffset, isCodingSegment, featureSpan } from '../model.js';

const NUM_W = 62;      // gutter for position numbers
const CHAR_W = 8.42;   // measured below; IBM Plex Mono at 13px
const GROUP = 10;

const LABEL_CH = 5.75; // width of one character at the 10px label size
const labelWidth = (text) => text.length * LABEL_CH + 8;

export function sequenceView(doc, opts = {}) {
  const { selection, onSelect, onPick, onPickSite, showTranslation = true, showSixFrame = false, showComplement = true, sites = [], highlight, reads = [], caret = null, overlays = [] } = opts;
  const L = doc.seq.length;
  // reading frames and primer sites ride along as extra bars
  const allFeatures = overlays.length ? [...doc.features, ...overlays] : doc.features;
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
    for (const f of allFeatures) {
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

  // Lay the cut labels for one line into rows, so that two sites a few bases
  // apart do not print over each other.
  function packSites(lineIdx) {
    const from = lineIdx * state.perLine;
    const to = Math.min(L, from + state.perLine);
    const ends = [];
    const placed = [];
    sites.forEach((st, idx) => {
      if (st.top <= from || st.top > to) return;
      const cx = NUM_W + (st.top - from) * charW;
      const w = labelWidth(st.enzyme.name);
      let row = 0;
      while (ends[row] !== undefined && ends[row] > cx) row++;
      ends[row] = cx + w;
      placed.push({ st, idx, cx, row });
    });
    return placed;
  }

  // a caret at the very end of a full last line belongs to that line, not a new one
  function caretLine() {
    if (caret === L && L > 0 && L % state.perLine === 0) return L / state.perLine - 1;
    return Math.floor(caret / state.perLine);
  }

  function buildLine(i) {
    const from = i * state.perLine;
    const to = Math.min(L, from + state.perLine);
    const top = doc.seq.slice(from, to);
    const line = h('div', { class: 'seqline', style: { position: 'absolute', left: 0, right: 0, top: `${i * state.lineH}px`, height: `${state.lineH}px` } });
    let y = 0;

    // enzyme cuts above the line, stacked so that close ones stay readable
    if (sites.length) {
      const rows = state.siteRows;
      const placed = packSites(i);
      const block = rows * 12 + 2;
      const strands = showComplement ? 34 : 17;
      for (const { st, idx, cx, row } of placed) {
        const top = (rows - 1 - row) * 12;
        line.append(h('div', {
          class: 'cutlabel', dataset: { site: String(idx) }, title: `${st.enzyme.name} cuts after ${st.top}`,
          style: { left: `${cx + 2}px`, top: `${top}px` },
        }, st.enzyme.name));
        line.append(h('div', {
          class: 'cutmark',
          style: { left: `${cx}px`, top: `${top + 11}px`, height: `${block + strands - top - 11}px` },
        }));
      }
      y += block;
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

    // the insertion point, between two bases
    if (caret !== null && caret >= 0 && caret <= L && caretLine() === i) {
      line.append(h('div', {
        class: 'caret',
        style: { left: `${NUM_W + (caret - from) * charW - 1}px`, top: `${y}px`, height: `${showComplement ? 34 : 17}px` },
      }));
    }

    line.append(h('div', { class: 'bases', style: { left: `${NUM_W}px`, top: `${y}px` } }, top));
    y += 17;
    if (showComplement) {
      line.append(h('div', { class: 'bases bottom', style: { left: `${NUM_W}px`, top: `${y}px` } }, complement(top)));
      y += 17;
    }

    // aligned sequencing reads, directly under the strands they align to
    for (const rd of reads) {
      const end = rd.start + rd.chars.length;
      if (end <= from || rd.start >= to) continue;
      const row = h('div', { class: 'readrow', style: { left: `${NUM_W}px`, top: `${y}px` } });
      const runs = [];
      let run = null;
      for (let p = from; p < to; p++) {
        let ch = ' ';
        let cls = 'blank';
        if (p >= rd.start && p < end) {
          const c = rd.chars[p - rd.start];
          if (c === '-') { ch = '-'; cls = 'del'; } else if (c === doc.seq[p]) { ch = '·'; cls = 'match'; } else { ch = c; cls = 'mis'; }
        }
        if (run && run.cls === cls) run.text += ch;
        else { run = { cls, text: ch }; runs.push(run); }
      }
      for (const r of runs) row.append(h('span', { class: r.cls }, r.text));
      line.append(h('div', { class: 'readname', style: { top: `${y}px` }, title: `${rd.label}${rd.reverse ? ' (reverse)' : ''}` },
        `${rd.reverse ? '◀' : '▶'} ${rd.label}`));
      line.append(row);
      for (const [at, bases] of rd.insertions) {
        if (at < from || at >= to) continue;
        line.append(h('div', {
          class: 'insmark', title: `${bases.length} base${bases.length > 1 ? 's' : ''} inserted: ${bases}`,
          style: { left: `${NUM_W + (at - from) * charW - 1}px`, top: `${y}px` },
        }));
      }
      y += 15;
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

    // six-frame translation track: top frames are codon-aligned; the three
    // reverse frames are shown as a compact block so they stay readable.
    if (showSixFrame) {
      for (let f = 1; f <= 3; f++) {
        const start = Math.max(0, from - mod(from - (f - 1), 3));
        const protein = translate(doc.seq.slice(start, Math.min(L, to + 3)));
        if (!protein) continue;
        line.append(h('div', {
          class: 'aa', style: { left: `${NUM_W + (start - from) * charW}px`, top: `${y}px`, opacity: .85 },
          title: `frame +${f}`,
        }, `+${f} ${[...protein].join('  ')}`));
        y += 13;
      }
      for (let f = 1; f <= 3; f++) {
        const rFrom = L - to;
        const start = Math.max(0, rFrom - mod(rFrom - (f - 1), 3));
        const stop = Math.min(L, rFrom + (to - from) + 3);
        // the reverse strand over [start, stop) is the reverse complement of the top strand over [L - stop, L - start);
        // doing the whole sequence here, once per line, was what made six-frame view slow on long sequences
        const protein = translate(revcomp(doc.seq.slice(L - stop, L - start)));
        if (!protein) continue;
        line.append(h('div', {
          class: 'aa', style: { left: `${NUM_W}px`, top: `${y}px`, opacity: .7 },
          title: `frame −${f}`,
        }, `−${f} ${[...protein].join('  ')}`));
        y += 13;
      }
    }

    // feature bars; a name too wide for its bar is printed beside it instead of
    // being clipped
    const lineWidth = NUM_W + state.perLine * charW;
    for (const { f, s0, s1, startsHere, gap } of featuresOnLine(from, to)) {
      const left = NUM_W + (s0 - from) * charW;
      const width = Math.max(2, (s1 - s0) * charW);
      const dim = highlight && highlight !== f.id ? 0.4 : 1;
      const bar = h('div', {
        class: 'fbar',
        title: f.name,
        dataset: { feature: f.id },
        style: {
          left: `${left}px`, width: `${width}px`, top: `${y}px`,
          background: gap || f.orf ? 'transparent' : f.color,
          borderTop: gap ? '2px dotted var(--muted)' : 'none',
          ...(f.orf ? { border: '1px dashed var(--muted)', color: 'var(--muted)', height: '13px' } : {}),
          opacity: dim,
        },
      });
      line.append(bar);
      if (!gap && startsHere) {
        const text = `${f.strand === -1 ? '◀ ' : ''}${f.name}${f.strand === 1 ? ' ▶' : ''}`;
        const need = labelWidth(text);
        if (need <= width) {
          bar.textContent = text;
        } else {
          const after = left + width + 4;
          const x = after + need <= lineWidth ? after : Math.max(NUM_W, left - need - 2);
          line.append(h('div', {
            class: 'foutside', title: f.name, dataset: { feature: f.id },
            style: { left: `${x}px`, top: `${y}px`, color: f.orf ? 'var(--muted)' : f.color, opacity: dim },
          }, text));
        }
      }
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
    for (const f of allFeatures) {
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
    let siteRows = 0;
    if (sites.length) {
      for (let i = 0; i < total; i++) siteRows = Math.max(siteRows, packSites(i).reduce((n, p) => Math.max(n, p.row + 1), 0));
    }
    state.siteRows = Math.max(1, siteRows);
    let maxReads = 0;
    for (let i = 0; i < total; i++) {
      const from = i * state.perLine;
      const to = from + state.perLine;
      let n = 0;
      for (const rd of reads) if (rd.start < to && rd.start + rd.chars.length > from) n++;
      maxReads = Math.max(maxReads, n);
    }
    const strands = 1 + (showComplement ? 1 : 0);
    return Math.round((sites.length ? state.siteRows * 12 + 2 : 0) + strands * 17 + maxReads * 15 + maxAA * 14 + (showSixFrame ? 6 * 13 : 0) + maxBars * 15 + 10);
  }

  function layout() {
    charW = probe.getBoundingClientRect().width / 100 || CHAR_W;
    const avail = wrap.clientWidth - NUM_W - 24;
    const perLine = Math.max(GROUP, Math.floor(avail / charW / GROUP) * GROUP);
    const was = state.perLine;
    state.perLine = perLine;           // lineHeight measures against this width
    const lineH = lineHeight();
    if (perLine === was && lineH === state.lineH && state.first >= 0) return false;
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

  let downOn = null;
  wrap.addEventListener('pointerdown', (ev) => {
    if (ev.button !== 0) return;
    downOn = ev.target.closest?.('[data-feature], [data-site]') ?? null;
    anchor = posAt(ev);
    // capture keeps the drag alive outside the element; not every pointer can be
    // captured, and failing to is not a reason to drop the selection
    try { wrap.setPointerCapture(ev.pointerId); } catch { /* carry on uncaptured */ }
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
    const clicked = downOn;
    anchor = null;
    downOn = null;
    // a click that did not drag, on a feature bar, its name or a cut label,
    // selects that thing rather than clearing the selection
    if (!next && clicked) {
      const fid = clicked.getAttribute('data-feature');
      if (fid) {
        const f = allFeatures.find((x) => x.id === fid);
        if (f && onPick) { onPick(f); return; }
      }
      const idx = clicked.getAttribute('data-site');
      if (idx !== null && sites[Number(idx)] && onPickSite) { onPickSite(sites[Number(idx)]); return; }
    }
    wrap.setSelection(next);
    // a click that did not drag leaves the insertion point there
    onSelect?.(next, { live: false, caret: next ? null : posAt(ev) });
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
