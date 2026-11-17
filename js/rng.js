// Small deterministic PRNG for reproducible simulations and property tests.
// mulberry32: 32-bit state, good enough for choosing among alternatives.

export function mulberry32(seed = 1) {
  let a = seed >>> 0;
  return function next() {
    a = (a + 0x6d2b79f5) >>> 0;
    let t = a;
    t = Math.imul(t ^ (t >>> 15), t | 1);
    t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

export function randomInt(rng, min, max) {
  return min + Math.floor(rng() * (max - min + 1));
}

export function pick(rng, list) {
  return list[Math.floor(rng() * list.length) % list.length];
}

export function randomBases(rng, length) {
  const bases = 'ACGT';
  let out = '';
  for (let i = 0; i < length; i++) out += bases[Math.floor(rng() * 4) % 4];
  return out;
}
