(function (root) {
  "use strict";
  const A = typeof module === "object" && module.exports ? require("../TeachingTools/statistics.js") : root.Astro;
  const defaults = { limit: 20, overlay: "show", seed: 42 };
  function observed(logLuminosity, distance, limit) {
    A.check([logLuminosity, distance, limit].every(Number.isFinite) && distance > 0,
      "Luminosity and magnitude must be finite; distance must be positive.");
    return -18 - 2.5 * logLuminosity + 5 * Math.log10(distance) + 25 <= limit;
  }
  function threshold(distance, limit) {
    A.check(Number.isFinite(distance) && distance > 0 && Number.isFinite(limit),
      "Distance must be positive and the magnitude limit finite.");
    return (5 * Math.log10(distance) + 7 - limit) / 2.5;
  }
  function atPhase(result, phase = 1) {
    A.check(Number.isFinite(phase) && phase >= 0 && phase <= 1, "Animation phase must lie between zero and one.");
    const limit = 26 + phase * (result.targetLimit - 26);
    const detected = result.distance.map((distance, i) => observed(result.luminosity[i], distance, limit));
    const sums = Array(10).fill(0), counts = Array(10).fill(0);
    let count = 0, sum = 0;
    detected.forEach((yes, i) => {
      if (!yes) return;
      count++; sum += result.luminosity[i];
      const bin = Math.min(9, Math.ceil(result.distance[i] / 80) - 1);
      sums[bin] += result.luminosity[i]; counts[bin]++;
    });
    return { limit, detected, count, missing: detected.length - count, completeness: count / detected.length,
      mean: count ? sum / count : null, binCounts: counts,
      binMeans: sums.map((total, i) => counts[i] ? total / counts[i] : null) };
  }
  function compute(parameters = {}) {
    const p = { ...defaults, ...parameters };
    A.check(Number.isFinite(p.limit) && p.limit >= 17 && p.limit <= 24, "Magnitude limit must be between 17 and 24.");
    A.check(["show", "hide"].includes(p.overlay), "Choose whether to show or hide nondetections.");
    const rng = A.createRng(p.seed), n = 2400;
    const distance = Array.from({ length: n }, () => Math.cbrt(30 ** 3 + rng.uniform() * (800 ** 3 - 30 ** 3)));
    const luminosity = Array.from({ length: n }, () => 0.6 * rng.normal());
    const result = { distance, luminosity, targetLimit: p.limit };
    return { ...result, final: atPhase(result) };
  }
  const api = { defaults, observed, threshold, atPhase, compute };
  if (typeof module === "object" && module.exports) module.exports = api;
  else root.SelectionEffectsModel = api;
})(globalThis);
