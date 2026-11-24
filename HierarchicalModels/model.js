(function (root) {
  "use strict";
  const A = typeof module === "object" && module.exports ? require("../TeachingTools/statistics.js") : root.Astro;
  const defaults = { pooling: "partial", scatter: 1, noise: 1, seed: 42 };
  function normalHierarchy(x, error, scatter, priorMean = 0, priorSd = 5) {
    A.check(Array.isArray(x) && x.length > 0 && x.every(Number.isFinite) &&
      Array.isArray(error) && error.length === x.length && error.every(s => Number.isFinite(s) && s > 0),
    "Supply equally sized finite observations and positive measurement errors.");
    A.check(Number.isFinite(scatter) && scatter >= 0 && Number.isFinite(priorMean) &&
      Number.isFinite(priorSd) && priorSd > 0, "Scatter must be nonnegative and the proper prior SD positive.");
    const totalVariance = error.map(s => scatter ** 2 + s ** 2);
    const populationVariance = 1 / (1 / priorSd ** 2 + totalVariance.reduce((sum, v) => sum + 1 / v, 0));
    const populationMean = populationVariance * (priorMean / priorSd ** 2 +
      x.reduce((sum, value, i) => sum + value / totalVariance[i], 0));
    const weight = totalVariance.map(v => scatter ** 2 / v);
    const mean = x.map((value, i) => weight[i] * value + (1 - weight[i]) * populationMean);
    const variance = error.map((s, i) => scatter ** 2 * s ** 2 / totalVariance[i] +
      (1 - weight[i]) ** 2 * populationVariance);
    A.check([populationMean, populationVariance, ...mean, ...variance].every(Number.isFinite) &&
      populationVariance > 0, "The hierarchy exceeds numerical precision.");
    return { mean, variance, populationMean, populationSd: Math.sqrt(populationVariance), weight };
  }
  function atPhase(result, phase = 1) {
    A.check(Number.isFinite(phase) && phase >= 0 && phase <= 1, "Animation phase must lie between zero and one.");
    return {
      estimate: result.observed.map((value, i) => value + phase * (result.posterior.mean[i] - value)),
      sd: result.error.map((s, i) => s + phase * (Math.sqrt(result.posterior.variance[i]) - s))
    };
  }
  function compute(parameters = {}) {
    const p = { ...defaults, ...parameters };
    A.check(["partial", "none", "complete"].includes(p.pooling), "Choose partial, no or complete pooling.");
    A.check(Number.isFinite(p.scatter) && p.scatter >= 0 && p.scatter <= 3, "Intrinsic scatter must be between 0 and 3 km/s.");
    A.check(Number.isFinite(p.noise) && p.noise >= 0.3 && p.noise <= 3, "Measurement multiplier must be between 0.3 and 3.");
    const rng = A.createRng(p.seed), n = 14;
    const truth = Array.from({ length: n }, () => 1.5 + p.scatter * rng.normal());
    const error = Array.from({ length: n }, (_, i) => (0.25 + 1.45 * i / (n - 1)) * p.noise);
    const observed = truth.map((value, i) => value + error[i] * rng.normal());
    const posterior = p.pooling === "none" ?
      { mean: [...observed], variance: error.map(s => s * s), weight: error.map(() => 1),
        populationMean: null, populationSd: null } :
      normalHierarchy(observed, error, p.pooling === "complete" ? 0 : p.scatter);
    const result = { truth, error, observed, posterior };
    return { ...result, final: atPhase(result) };
  }
  const api = { defaults, normalHierarchy, atPhase, compute };
  if (typeof module === "object" && module.exports) module.exports = api;
  else root.HierarchicalModelsModel = api;
})(globalThis);
