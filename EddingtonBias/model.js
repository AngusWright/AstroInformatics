(function (root) {
  "use strict";
  const A = typeof module === "object" && module.exports ? require("../TeachingTools/statistics.js") : root.Astro;
  const defaults = { noise: 0.3, slope: 1.5, threshold: 10.5, seed: 42 };
  function massDensity(x, slope, noise = 0, lower = 8.5, upper = 12) {
    A.check([x, slope, noise, lower, upper].every(Number.isFinite) && slope > 0 && noise >= 0 && upper > lower,
      "The density needs a positive slope, nonnegative noise and finite ordered bounds.");
    const normalisation = -Math.expm1(-slope * (upper - lower));
    if (noise === 0) return x >= lower && x <= upper ? slope * Math.exp(-slope * (x - lower)) / normalisation : 0;
    const a = (lower - x + slope * noise * noise) / noise;
    const b = (upper - x + slope * noise * noise) / noise;
    // Reflect positive tails to avoid subtracting two CDFs rounded to one.
    const interval = a > 0 ? A.normalCdf(-a) - A.normalCdf(-b) : A.normalCdf(b) - A.normalCdf(a);
    const value = interval > 0 ? Math.exp(Math.log(slope / normalisation) +
      slope * (lower - x) + (slope * noise) ** 2 / 2 + Math.log(interval)) : 0;
    A.check(Number.isFinite(value), "The density exceeds numerical precision.");
    return value;
  }
  function atPhase(result, phase = 1) {
    A.check(Number.isFinite(phase) && phase >= 0 && phase <= 1, "Animation phase must lie between zero and one.");
    const noise = phase * result.targetNoise;
    const observed = result.truth.map((value, i) => value + noise * result.errors[i]);
    let up = 0, down = 0, trueHigh = 0, observedHigh = 0;
    observed.forEach((value, i) => {
      const before = result.truth[i] >= result.threshold, after = value >= result.threshold;
      if (before) trueHigh++;
      if (after) observedHigh++;
      if (!before && after) up++;
      if (before && !after) down++;
    });
    return { observed, noise, up, down, net: up - down, trueHigh, observedHigh,
      lowLow: observed.length - trueHigh - up, highHigh: trueHigh - down,
      histogram: A.histogram(observed, result.breaks) };
  }
  function compute(parameters = {}) {
    const p = { ...defaults, ...parameters };
    A.check(Number.isFinite(p.noise) && p.noise >= 0 && p.noise <= 0.6, "Measurement error must be between 0 and 0.6 dex.");
    A.check(Number.isFinite(p.slope) && p.slope >= 0.3 && p.slope <= 3, "Population slope must be between 0.3 and 3.");
    A.check(Number.isFinite(p.threshold) && p.threshold >= 9.5 && p.threshold <= 11.5,
      "The high-mass threshold must be between 9.5 and 11.5 dex.");
    const rng = A.createRng(p.seed), n = 8000;
    const truth = Array.from({ length: n }, () => 8.5 - Math.log1p(-rng.uniform() * -Math.expm1(-p.slope * 3.5)) / p.slope);
    const errors = Array.from({ length: n }, () => rng.normal());
    const finalObserved = truth.map((value, i) => value + p.noise * errors[i]);
    const lo = Math.floor(Math.min(8.5, ...finalObserved) * 10) / 10 - 0.1;
    const hi = Math.ceil(Math.max(12, ...finalObserved) * 10) / 10 + 0.1;
    const bins = Math.round((hi - lo) * 10);
    const breaks = Array.from({ length: bins + 1 }, (_, i) => lo + i / 10);
    const result = { truth, errors, breaks, targetNoise: p.noise, threshold: p.threshold,
      slope: p.slope, trueHistogram: A.histogram(truth, breaks) };
    return { ...result, final: atPhase(result) };
  }
  const api = { defaults, massDensity, atPhase, compute };
  if (typeof module === "object" && module.exports) module.exports = api;
  else root.EddingtonBiasModel = api;
})(globalThis);
