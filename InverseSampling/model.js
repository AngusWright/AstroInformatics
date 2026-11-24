(function (root, factory) {
  "use strict";
  const model = factory(typeof module === "object" && module.exports ? require("../TeachingTools/statistics.js") : root.Astro);
  if (typeof module === "object" && module.exports) module.exports = model;
  else root.InverseSamplingModel = model;
})(globalThis, function (A) {
  "use strict";
  function probability(p) {
    A.check(Number.isFinite(p) && p >= 0 && p <= 1, "A probability must be between zero and one.");
  }
  function distribution(family, slope = 2.35) {
    A.check(["power", "exponential", "lognormal"].includes(family), "Choose a supported positive target distribution.");
    A.check(Number.isFinite(slope) && slope >= 0.5 && slope <= 3.5, "Power-law exponent must be between 0.5 and 3.5.");
    if (family === "exponential") return {
      label: "Exponential, rate 1",
      cdf: x => x <= 0 ? 0 : -Math.expm1(-x),
      density: x => x < 0 ? 0 : Math.exp(-x),
      quantile(p) { probability(p); return -Math.log1p(-p); }
    };
    if (family === "lognormal") return {
      label: "Log-normal, log-mean 0, log-SD 0.7",
      cdf: x => x <= 0 ? 0 : x === Infinity ? 1 : A.normalCdf(Math.log(x) / 0.7),
      density: x => x <= 0 || x === Infinity ? 0 : A.normalPdf(Math.log(x) / 0.7) / (0.7 * x),
      quantile(p) { probability(p); return Math.exp(0.7 * A.normalQuantile(p)); }
    };
    const lower = 0.1, upper = 20, logRange = Math.log(upper / lower), k = 1 - slope;
    const denominator = Math.expm1(k * logRange);
    const normalisation = k === 0 ? logRange : Math.exp(k * Math.log(lower)) * denominator / k;
    return {
      label: `Power law, alpha = ${slope}`, lower, upper,
      cdf(x) {
        if (x <= lower) return 0;
        if (x >= upper) return 1;
        return k === 0 ? Math.log(x / lower) / logRange : Math.expm1(k * Math.log(x / lower)) / denominator;
      },
      density: x => x < lower || x > upper ? 0 : Math.pow(x, -slope) / normalisation,
      quantile(p) {
        probability(p);
        if (p === 0) return lower;
        if (p === 1) return upper;
        return lower * Math.exp(k === 0 ? p * logRange : Math.log1p(p * denominator) / k);
      }
    };
  }
  function validate(p) {
    distribution(p.family, p.slope);
    A.check(Number.isInteger(p.seed) && p.seed >= 0 && p.seed <= 2147483647, "Seed must be an integer from 0 to 2147483647.");
    A.check(Number.isInteger(p.n) && p.n >= 50 && p.n <= 1500, "Use between 50 and 1500 random draws.");
  }
  function prepare(seed) {
    A.check(Number.isInteger(seed) && seed >= 0 && seed <= 2147483647, "Seed must be an integer from 0 to 2147483647.");
    const rng = A.createRng(seed);
    return { seed, uniforms: Array.from({ length: 1500 }, () => rng.uniform()) };
  }
  function ksStatistic(probabilities) {
    A.check(probabilities.every(p => Number.isFinite(p) && p >= 0 && p <= 1), "CDF values must be valid probabilities.");
    if (!probabilities.length) return null;
    const ordered = [...probabilities].sort((a, b) => a - b), n = ordered.length;
    return Math.max(...ordered.map((p, i) => Math.max((i + 1) / n - p, p - i / n)));
  }
  function compute(parameters, prepared) {
    validate(parameters);
    const data = prepared || prepare(parameters.seed), target = distribution(parameters.family, parameters.slope);
    A.check(data.seed === parameters.seed && data.uniforms.length >= parameters.n, "Cached uniform draws do not match these settings.");
    const uniforms = data.uniforms.slice(0, parameters.n), values = uniforms.map(target.quantile);
    const limits = [Math.log10(target.quantile(0.00001)), Math.log10(target.quantile(0.99999))];
    values.forEach(x => { limits[0] = Math.min(limits[0], Math.log10(x)); limits[1] = Math.max(limits[1], Math.log10(x)); });
    if (parameters.family === "power") { limits[0] = Math.log10(0.1); limits[1] = Math.log10(20); }
    return { data, parameters: { ...parameters }, target, uniforms, values, limits };
  }
  function atFrame(result, frame = 100) {
    A.check(Number.isInteger(frame) && frame >= 0 && frame <= 100, "Animation frame must be from 0 to 100.");
    const count = Math.floor(frame / 100 * result.parameters.n);
    const uniforms = result.uniforms.slice(0, count), values = result.values.slice(0, count);
    return { count, uniforms, values, discrepancy: ksStatistic(uniforms) };
  }
  function logDensity(target, z) {
    const x = z === Math.log10(target.lower) ? target.lower :
      z === Math.log10(target.upper) ? target.upper : Math.pow(10, z);
    return Math.LN10 * x * target.density(x);
  }
  return { distribution, validate, prepare, ksStatistic, compute, atFrame, logDensity };
});
