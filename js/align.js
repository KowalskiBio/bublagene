// Pairwise alignment with affine gaps: global, local and banded.
//
// Defaults: match +2, mismatch −1, gap open 3, gap extend 1 (a gap of length k
// costs open + k·extend). Results are ordinary objects so they can be recorded,
// compared and printed.

const DEFAULTS = { match: 2, mismatch: -1, gapOpen: 3, gapExtend: 1 };

function normalise(seq) {
  return String(seq).toUpperCase();
}

// Full DP with three states. `band` restricts |i-j| <= band (useful for long,
// similar sequences); pass null for unbanded.
function dp(a, b, { match, mismatch, gapOpen, gapExtend, band = null, local = false }) {
  const n = a.length;
  const m = b.length;
  const NEG = -1e15;
  const w = m + 1;
  const M = new Float64Array((n + 1) * w).fill(NEG);
  const X = new Float64Array((n + 1) * w).fill(NEG); // gap in b (consume a)
  const Y = new Float64Array((n + 1) * w).fill(NEG); // gap in a (consume b)
  const ptrM = new Int8Array((n + 1) * w);
  const ptrX = new Int8Array((n + 1) * w);
  const ptrY = new Int8Array((n + 1) * w);
  const idx = (i, j) => i * w + j;
  const allowed = (i, j) => band == null || Math.abs(i - j) <= band;
  M[idx(0, 0)] = 0;
  for (let i = 1; i <= n; i++) {
    if (!allowed(i, 0)) continue;
    M[idx(i, 0)] = local ? 0 : NEG;
    X[idx(i, 0)] = -(gapOpen + i * gapExtend);
    ptrX[idx(i, 0)] = i > 1 ? 1 : 0;
  }
  for (let j = 1; j <= m; j++) {
    if (!allowed(0, j)) continue;
    M[idx(0, j)] = local ? 0 : NEG;
    Y[idx(0, j)] = -(gapOpen + j * gapExtend);
    ptrY[idx(0, j)] = j > 1 ? 1 : 0;
  }
  let best = { score: local ? 0 : NEG, i: 0, j: 0, state: local ? 'M' : null };
  for (let i = 1; i <= n; i++) {
    const jStart = band == null ? 1 : Math.max(1, i - band);
    const jEnd = band == null ? m : Math.min(m, i + band);
    for (let j = jStart; j <= jEnd; j++) {
      const k = idx(i, j);
      const s = a[i - 1] === b[j - 1] ? match : mismatch;
      const diag = [M[idx(i - 1, j - 1)], X[idx(i - 1, j - 1)], Y[idx(i - 1, j - 1)]];
      const dmax = Math.max(...diag);
      let mv = dmax + s;
      ptrM[k] = diag[0] >= diag[1] && diag[0] >= diag[2] ? 0 : diag[1] >= diag[2] ? 1 : 2;
      if (local && mv < 0) { mv = 0; ptrM[k] = 3; }
      M[k] = mv;
      const openFromM = M[idx(i - 1, j)] - gapOpen - gapExtend;
      const openFromX = X[idx(i - 1, j)] - gapExtend;
      if (openFromM >= openFromX) { X[k] = openFromM; ptrX[k] = 0; } else { X[k] = openFromX; ptrX[k] = 1; }
      const openFromM2 = M[idx(i, j - 1)] - gapOpen - gapExtend;
      const openFromY = Y[idx(i, j - 1)] - gapExtend;
      if (openFromM2 >= openFromY) { Y[k] = openFromM2; ptrY[k] = 0; } else { Y[k] = openFromY; ptrY[k] = 1; }
      if (local) {
        if (X[k] < 0) X[k] = NEG;
        if (Y[k] < 0) Y[k] = NEG;
      }
      if (M[k] > best.score) best = { score: M[k], i, j, state: 'M' };
      if (X[k] > best.score) best = { score: X[k], i, j, state: 'X' };
      if (Y[k] > best.score) best = { score: Y[k], i, j, state: 'Y' };
    }
  }
  if (!local) {
    const end = Math.max(M[idx(n, m)], X[idx(n, m)], Y[idx(n, m)]);
    const state = end === M[idx(n, m)] ? 'M' : end === X[idx(n, m)] ? 'X' : 'Y';
    best = { score: end, i: n, j: m, state };
  }
  return { a, b, M, X, Y, ptrM, ptrX, ptrY, w, best, local };
}

function traceback(d) {
  const { a, b, w, ptrM, ptrX, ptrY, best } = d;
  const idx = (i, j) => i * w + j;
  let i = best.i;
  let j = best.j;
  let state = best.state;
  const outA = [];
  const outB = [];
  const ops = [];
  while (i > 0 || j > 0) {
    if (d.local && state === 'M' && d.M[idx(i, j)] === 0) break;
    if (state === 'M') {
      outA.push(a[i - 1]);
      outB.push(b[j - 1]);
      ops.push(a[i - 1] === b[j - 1] ? '=' : 'X');
      const p = ptrM[idx(i, j)];
      if (p === 3) break;
      i--; j--;
      state = p === 0 ? 'M' : p === 1 ? 'X' : 'Y';
    } else if (state === 'X') {
      outA.push(a[i - 1]);
      outB.push('-');
      ops.push('D');
      const p = ptrX[idx(i, j)];
      i--;
      state = p === 0 ? 'M' : 'X';
    } else {
      outA.push('-');
      outB.push(b[j - 1]);
      ops.push('I');
      const p = ptrY[idx(i, j)];
      j--;
      state = p === 0 ? 'M' : 'Y';
    }
  }
  return {
    alignedA: outA.reverse().join(''),
    alignedB: outB.reverse().join(''),
    ops: ops.reverse().join(''),
    startA: best.i - outA.length + outA.filter((c) => c === '-').length,
    endA: best.i,
    startB: best.j - outB.length + outB.filter((c) => c === '-').length,
    endB: best.j,
  };
}

export function alignPair(seqA, seqB, opts = {}) {
  const a = normalise(seqA);
  const b = normalise(seqB);
  const options = { ...DEFAULTS, band: null, local: false, ...opts };
  const d = dp(a, b, options);
  const t = traceback(d);
  let matches = 0;
  let mismatches = 0;
  let gaps = 0;
  let columns = 0;
  for (let i = 0; i < t.alignedA.length; i++) {
    columns++;
    if (t.alignedA[i] === '-' || t.alignedB[i] === '-') gaps++;
    else if (t.alignedA[i] === t.alignedB[i]) matches++;
    else mismatches++;
  }
  return {
    ...t,
    score: d.best.score,
    matches,
    mismatches,
    gaps,
    columns,
    identity: columns ? (100 * matches) / columns : 0,
    identityAligned: matches + mismatches ? (100 * matches) / (matches + mismatches) : 0,
    global: !options.local,
    band: options.band,
    options: { match: options.match, mismatch: options.mismatch, gapOpen: options.gapOpen, gapExtend: options.gapExtend },
  };
}

export function alignGlobal(a, b, opts = {}) {
  return alignPair(a, b, { ...opts, local: false });
}

export function alignLocal(a, b, opts = {}) {
  return alignPair(a, b, { ...opts, local: true });
}

export function alignBanded(a, b, { band = 8, ...opts } = {}) {
  return alignPair(a, b, { ...opts, local: false, band });
}

// Dot-plot points: one per window pair whose identity clears `threshold`.
export function dotPlot(a, b, { window = 10, step = 1, threshold = 0.7, maxPoints = 20000 } = {}) {
  const x = normalise(a);
  const y = normalise(b);
  const points = [];
  let best = 0;
  for (let i = 0; i + window <= x.length; i += step) {
    for (let j = 0; j + window <= y.length; j += step) {
      let same = 0;
      for (let k = 0; k < window; k++) if (x[i + k] === y[j + k]) same++;
      const score = same / window;
      if (score > best) best = score;
      if (score >= threshold) {
        points.push([i, j, Math.round(score * 1000) / 1000]);
        if (points.length >= maxPoints) return { points, best, width: x.length, height: y.length, truncated: true, window };
      }
    }
  }
  return { points, best, width: x.length, height: y.length, truncated: false, window };
}

// Character grid of the dots (handy for the CLI).
export function dotPlotText(a, b, { window = 6, step = 2, threshold = 0.8, width = 60 } = {}) {
  const plot = dotPlot(a, b, { window, step, threshold });
  const cols = Math.min(width, Math.ceil(plot.width / step));
  const rows = Math.min(width, Math.ceil(plot.height / step));
  const grid = Array.from({ length: rows }, () => Array(cols).fill(' '));
  for (const [i, j] of plot.points) {
    const c = Math.floor(i / step);
    const r = Math.floor(j / step);
    if (r < rows && c < cols) grid[r][c] = '•';
  }
  return grid.map((row) => row.join('')).join('\n');
}

export function formatAlignment(result, { width = 60 } = {}) {
  const lines = [];
  const { alignedA, alignedB, ops } = result;
  const useOps = ops.length === alignedA.length;
  for (let i = 0; i < alignedA.length; i += width) {
    const a = alignedA.slice(i, i + width);
    const b = alignedB.slice(i, i + width);
    const mid = useOps
      ? ops.slice(i, i + width).replace(/=/g, '|').replace(/X/g, ' ').replace(/[DI]/g, ' ')
      : [...a].map((c, k) => (c === b[k] && c !== '-' ? '|' : ' ')).join('');
    lines.push(`  ${a}`);
    lines.push(`  ${mid}`);
    lines.push(`  ${b}`);
    lines.push('');
  }
  lines.push(`${result.matches} matches, ${result.mismatches} mismatches, ${result.gaps} gap columns, identity ${result.identity.toFixed(1)}%, score ${result.score}`);
  return lines.join('\n');
}
