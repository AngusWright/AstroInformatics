(function (root) {
  "use strict";
  const A = typeof module === "object" && module.exports ? require("../TeachingTools/statistics.js") : root.Astro;
  function predictiveNormal(x, intrinsic, measurement, priorSd = 5) {
    A.check(Array.isArray(x) && x.length > 0 && x.every(Number.isFinite), "Provide at least one finite observed velocity.");
    A.check(Number.isFinite(intrinsic) && intrinsic >= 0, "Intrinsic dispersion must be non-negative.");
    A.check(Number.isFinite(measurement) && measurement > 0, "Measurement SD must be positive.");
    A.check(Number.isFinite(priorSd) && priorSd > 0, "Prior SD must be positive.");
    const observationVariance = intrinsic ** 2 + measurement ** 2;
    const meanVariance = 1 / (1 / priorSd ** 2 + x.length / observationVariance);
    const mean = meanVariance * x.reduce((sum, value) => sum + value, 0) / observationVariance;
    const variance = { mean: meanVariance, star: meanVariance + intrinsic ** 2,
      measurement: meanVariance + observationVariance };
    A.check(Number.isFinite(mean) && Object.values(variance).every(v => Number.isFinite(v) && v > 0),
      "The requested model exceeds numeric precision.");
    const z = A.normalQuantile(0.975);
    const intervals = Object.fromEntries(Object.entries(variance).map(([name, value]) =>
      [name, [mean - z * Math.sqrt(value), mean + z * Math.sqrt(value)]]));
    return { mean, variance, intervals, observationVariance };
  }
  function compute(p) {
    A.check(Number.isInteger(p.seed) && p.seed >= 0 && p.seed <= 2147483647, "Seed must be an integer from 0 to 2147483647.");
    A.check(Number.isInteger(p.n) && p.n >= 1 && p.n <= 300, "Sample size must be an integer from 1 to 300.");
    A.check(Number.isFinite(p.intrinsic) && p.intrinsic >= 0 && p.intrinsic <= 4, "Intrinsic dispersion must be between 0 and 4.");
    A.check(Number.isFinite(p.measurement) && p.measurement >= 0.1 && p.measurement <= 3, "Measurement SD must be between 0.1 and 3.");
    A.check(typeof p.show_data === "boolean", "Observation display must be on or off.");
    const rng = A.createRng(p.seed);
    const stars = Array.from({ length: 300 }, () => 2 + p.intrinsic * rng.normal());
    const observations = stars.map(star => star + p.measurement * rng.normal());
    const x = observations.slice(0, p.n);
    return { stars, observations, x, posterior: predictiveNormal(x, p.intrinsic, p.measurement) };
  }
  const api = { predictiveNormal, compute };
  if (typeof module === "object" && module.exports) module.exports = api;
  else root.PosteriorPredictiveModel = api;
})(globalThis);
