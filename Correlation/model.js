(function (root, factory) {
  "use strict";
  const api = factory(typeof module === "object" && module.exports ? require("../TeachingTools/statistics.js") : root.Astro);
  if (typeof module === "object" && module.exports) module.exports = api;
  root.CorrelationModel = api;
})(globalThis, function (A) {
  "use strict";
  const defaults = { seed: 42, scatter: 0.4, centre: 0.5, width: 0.2 }, SIZE = 1800;
  function validate(p) {
    A.check(Number.isInteger(p.seed) && p.seed >= 0 && p.seed <= 2147483647, "Enter an integer seed between 0 and 2147483647.");
    for (const [key, lo, hi] of [["scatter", 0, 2], ["centre", 0, 1], ["width", 0.05, 1]]) {
      A.check(Number.isFinite(p[key]) && p[key] >= lo && p[key] <= hi, `${key} must be between ${lo} and ${hi}.`);
    }
  }
  function prepare(seed) {
    const rng = A.createRng(seed);
    return { seed, distance: Array.from({ length: SIZE }, () => rng.uniform()),
      normalX: Array.from({ length: SIZE }, () => rng.normal()), normalY: Array.from({ length: SIZE }, () => rng.normal()) };
  }
  function statistics(x, y) {
    A.check(x.length === y.length && x.every(Number.isFinite) && y.every(Number.isFinite), "Correlation needs matching finite observations.");
    const n = x.length;
    if (!n) return { n, meanX: null, meanY: null, correlation: null, slope: null, intercept: null };
    const meanX = x[0] + x.reduce((sum, value) => sum + (value - x[0]) / n, 0);
    const meanY = y[0] + y.reduce((sum, value) => sum + (value - y[0]) / n, 0);
    let xx = 0, yy = 0, xy = 0;
    x.forEach((value, i) => { const dx = value - meanX, dy = y[i] - meanY; xx += dx * dx; yy += dy * dy; xy += dx * dy; });
    const denominator = Math.sqrt(xx) * Math.sqrt(yy);
    const correlation = n >= 3 && denominator > 0 ? Math.max(-1, Math.min(1, xy / denominator)) : null;
    const slope = n >= 3 && xx > 0 ? xy / xx : null;
    return { n, meanX, meanY, correlation, slope, intercept: slope === null ? null : meanY - slope * meanX };
  }
  function populationCorrelation(scatter, distanceWidth = 1) {
    A.check(Number.isFinite(scatter) && scatter >= 0 && Number.isFinite(distanceWidth) && distanceWidth >= 0 && distanceWidth <= 1,
      "Scatter must be non-negative and distance width between 0 and 1.");
    const commonVariance = distanceWidth ** 2 / 3, totalVariance = commonVariance + scatter ** 2;
    return totalVariance > 0 ? commonVariance / totalVariance : null;
  }
  function compute(parameters, prepared) {
    const p = { ...defaults, ...parameters };
    validate(p);
    const base = prepared || prepare(p.seed);
    A.check(base.seed === p.seed && base.distance.length === SIZE, "The cached population must match the seed.");
    const intrinsicX = base.normalX.map(z => p.scatter * z), intrinsicY = base.normalY.map(z => p.scatter * z);
    const x = base.distance.map((z, i) => -2 * z + intrinsicX[i]);
    const y = base.distance.map((z, i) => -2 * z + intrinsicY[i]);
    const selected = base.distance.map(z => Math.abs(z - p.centre) <= p.width / 2);
    const sliceX = x.filter((_, i) => selected[i]), sliceY = y.filter((_, i) => selected[i]);
    const lower = Math.max(0, p.centre - p.width / 2), upper = Math.min(1, p.centre + p.width / 2);
    return { base, distance: base.distance, x, y, intrinsicX, intrinsicY, selected,
      full: statistics(x, y), slice: statistics(sliceX, sliceY), intrinsic: statistics(intrinsicX, intrinsicY),
      theoryFull: populationCorrelation(p.scatter), theorySlice: populationCorrelation(p.scatter, upper - lower),
      lower, upper };
  }
  return { defaults, SIZE, validate, prepare, statistics, populationCorrelation, compute };
});
