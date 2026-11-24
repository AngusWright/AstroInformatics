(function (root, factory) {
  "use strict";
  const model = factory(typeof module === "object" && module.exports ? require("../TeachingTools/statistics.js") : root.Astro);
  if (typeof module === "object" && module.exports) module.exports = model;
  else root.DensityEstimationModel = model;
})(globalThis, function (A) {
  "use strict";
  function populationDensity(x) {
    return 0.6 * A.normalPdf(x, -0.8, 0.35) + 0.4 * A.normalPdf(x, 0.9, 0.35);
  }
  function validate(p) {
    A.check(Number.isInteger(p.seed) && p.seed >= 0 && p.seed <= 2147483647, "Seed must be an integer from 0 to 2147483647.");
    A.check(Number.isInteger(p.n) && p.n >= 50 && p.n <= 1000, "Use between 50 and 1000 observed objects.");
    A.check(Number.isFinite(p.width) && p.width >= 0.1 && p.width <= 1.5, "Bin width must be between 0.1 and 1.5.");
    A.check(Number.isFinite(p.origin) && p.origin >= 0 && p.origin <= 1, "Bin origin must be a fraction from zero to one.");
    A.check(Number.isFinite(p.bandwidth) && p.bandwidth >= 0.05 && p.bandwidth <= 1, "KDE bandwidth must be between 0.05 and 1.");
  }
  function prepare(seed) {
    A.check(Number.isInteger(seed) && seed >= 0 && seed <= 2147483647, "Seed must be an integer from 0 to 2147483647.");
    const rng = A.createRng(seed), components = Array.from({ length: 1000 }, () => rng.uniform() < 0.6);
    return { seed, values: components.map(first => (first ? -0.8 : 0.9) + 0.35 * rng.normal()) };
  }
  function checkKde(values, bandwidth) {
    A.check(values.length > 0 && values.every(Number.isFinite), "A KDE needs finite observations.");
    A.check(Number.isFinite(bandwidth) && bandwidth > 0, "KDE bandwidth must be positive.");
  }
  function kdeAt(values, bandwidth, x) {
    checkKde(values, bandwidth);
    A.check(Number.isFinite(x), "KDE evaluation position must be finite.");
    return values.reduce((sum, value) => sum + A.normalPdf(x, value, bandwidth), 0) / values.length;
  }
  function kdeMass(values, bandwidth, lower, upper) {
    checkKde(values, bandwidth);
    A.check(Number.isFinite(lower) && Number.isFinite(upper) && upper >= lower, "KDE mass needs ordered finite bounds.");
    return values.reduce((sum, value) => sum + A.normalCdf(upper, value, bandwidth) - A.normalCdf(lower, value, bandwidth), 0) / values.length;
  }
  function kdeOutsideMass(values, bandwidth, lower, upper) {
    checkKde(values, bandwidth);
    A.check(Number.isFinite(lower) && Number.isFinite(upper) && upper >= lower, "KDE tails need ordered finite bounds.");
    return values.reduce((sum, value) =>
      sum + A.normalCdf(lower, value, bandwidth) + A.normalCdf(-upper, -value, bandwidth), 0) / values.length;
  }
  function countPeaks(densities) {
    A.check(densities.length >= 3 && densities.every(y => Number.isFinite(y) && y >= 0), "Peak counting needs a nonnegative density grid.");
    const threshold = 0.02 * Math.max(...densities), peaks = [];
    for (let i = 1; i < densities.length - 1; i++) {
      if (densities[i] <= densities[i - 1]) continue;
      let end = i;
      while (end + 1 < densities.length && densities[end + 1] === densities[i]) end++;
      if (end < densities.length - 1 && densities[end] > densities[end + 1] && densities[i] > threshold) {
        peaks.push(Math.floor((i + end) / 2));
      }
      i = end;
    }
    return peaks;
  }
  function compute(parameters, prepared) {
    validate(parameters);
    const data = prepared || prepare(parameters.seed);
    A.check(data.seed === parameters.seed && data.values.length >= parameters.n, "Cached observations do not match these settings.");
    const values = data.values.slice(0, parameters.n), offset = parameters.origin * parameters.width;
    const minimum = Math.min(...values), maximum = Math.max(...values);
    let first = Math.floor((minimum - offset) / parameters.width);
    let last = Math.floor((maximum - offset) / parameters.width) + 1;
    if (first * parameters.width + offset > minimum) first--;
    if (last * parameters.width + offset < maximum) last++;
    const breaks = Array.from({ length: last - first + 1 }, (_, i) => (first + i) * parameters.width + offset);
    const histogram = A.histogram(values, breaks);
    A.check(histogram.excluded === 0, "The histogram bounds must include every observation.");
    const limits = [Math.min(-3, breaks[0]), Math.max(3, breaks[breaks.length - 1])];
    const x = Array.from({ length: 801 }, (_, i) => limits[0] + (limits[1] - limits[0]) * i / 800);
    const density = x.map(location =>
      values.reduce((sum, value) => sum + A.normalPdf(location, value, parameters.bandwidth), 0) / values.length);
    return { data, parameters: { ...parameters }, values, histogram, breaks, limits, x, density, peaks: countPeaks(density),
      histogramArea: histogram.density.reduce((sum, height, i) => sum + height * histogram.widths[i], 0),
      outsideMass: kdeOutsideMass(values, parameters.bandwidth, ...limits) };
  }
  return { populationDensity, validate, prepare, kdeAt, kdeMass, kdeOutsideMass, countPeaks, compute };
});
