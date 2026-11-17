// Plots along a sequence: GC content and GC skew in a sliding window.
// Circular sequences wrap, so a window near the origin is not cut short.

const at = (seq, i, circular) => {
  const L = seq.length;
  if (circular) return seq[((i % L) + L) % L];
  return i >= 0 && i < L ? seq[i] : '';
};

function profile(seq, window, step, circular, value) {
  const L = seq.length;
  if (!L) return [];
  const w = Math.max(2, Math.min(window, L));
  const out = [];
  const half = Math.floor(w / 2);
  for (let centre = 0; centre < L; centre += step) {
    let g = 0; let c = 0; let n = 0;
    for (let k = -half; k < w - half; k++) {
      const b = at(seq, centre + k, circular);
      if (!b) continue;
      n++;
      if (b === 'G') g++; else if (b === 'C') c++;
    }
    if (n) out.push({ position: centre, value: value(g, c, n) });
  }
  return out;
}

// Percent G+C in each window.
export const gcProfile = (seq, window = 100, step = 10, circular = false) =>
  profile(seq.toUpperCase(), window, step, circular, (g, c, n) => (100 * (g + c)) / n);

// (G − C) / (G + C) in each window: it changes sign at replication origins and termini in bacteria.
export const gcSkewProfile = (seq, window = 100, step = 10, circular = false) =>
  profile(seq.toUpperCase(), window, step, circular, (g, c) => (g + c ? (g - c) / (g + c) : 0));

// Running total of GC skew: its minima and maxima mark where the skew flips.
export function cumulativeSkew(seq, step = 10) {
  const up = seq.toUpperCase();
  const out = [];
  let total = 0;
  for (let i = 0; i < up.length; i++) {
    if (up[i] === 'G') total++; else if (up[i] === 'C') total--;
    if (i % step === 0) out.push({ position: i, value: total });
  }
  return out;
}

export function summarize(points) {
  if (!points.length) return { min: 0, max: 0, mean: 0 };
  let min = Infinity; let max = -Infinity; let sum = 0;
  for (const p of points) { if (p.value < min) min = p.value; if (p.value > max) max = p.value; sum += p.value; }
  return { min, max, mean: sum / points.length };
}
