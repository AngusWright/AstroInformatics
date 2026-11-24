(function (root, factory) {
  "use strict";
  const model = factory(typeof module === "object" && module.exports ? require("../TeachingTools/statistics.js") : root.Astro);
  if (typeof module === "object" && module.exports) module.exports = model;
  else root.RobustStatisticsModel = model;
})(globalThis, function (A) {
  "use strict";
  function trimmedMean(values, trim = 0.1) {
    A.check(values.length > 0 && values.every(Number.isFinite), "A trimmed mean needs finite observations.");
    A.check(Number.isFinite(trim) && trim >= 0 && trim < 0.5, "Trim a fraction from zero up to, but not including, one half per tail.");
    const ordered = [...values].sort((a, b) => a - b), cut = Math.floor(trim * ordered.length);
    const kept = ordered.slice(cut, ordered.length - cut);
    return kept.reduce((sum, value) => sum + value, 0) / kept.length;
  }
  function summaries(values, trim = 0.1) {
    A.check(values.length >= 2 && values.every(Number.isFinite), "Summaries need at least two finite measurements.");
    const sample = A.summarize(values), median = A.quantile(values, 0.5);
    return { mean: sample.mean, median, sd: sample.sd,
      nmad: A.quantile(values.map(value => Math.abs(value - median)), 0.5) / A.normalQuantile(0.75),
      trimmed: trimmedMean(values, trim) };
  }
  function validate(p) {
    A.check(Number.isInteger(p.seed) && p.seed >= 0 && p.seed <= 2147483647, "Seed must be an integer from 0 to 2147483647.");
    A.check(Number.isInteger(p.n) && p.n >= 20 && p.n <= 120, "Use between 20 and 120 measurements.");
    A.check(Number.isInteger(p.contaminants) && p.contaminants >= 0 && p.contaminants <= 10,
      "Replace between zero and ten measurements.");
    A.check(Number.isFinite(p.position) && p.position >= -12 && p.position <= 12, "Contaminant velocity must be between -12 and 12 km/s.");
    A.check(Number.isFinite(p.trim) && p.trim >= 0 && p.trim <= 0.4, "Trim between zero and 0.4 of each tail.");
  }
  function prepare(seed) {
    A.check(Number.isInteger(seed) && seed >= 0 && seed <= 2147483647, "Seed must be an integer from 0 to 2147483647.");
    const rng = A.createRng(seed);
    return { seed, clean: Array.from({ length: 120 }, () => rng.normal()),
      offsets: Array.from({ length: 10 }, () => 0.1 * rng.normal()),
      jitter: Array.from({ length: 120 }, () => 2.25 + 0.45 * rng.uniform()) };
  }
  function contaminate(clean, offsets, count, position) {
    A.check(Number.isInteger(count) && count >= 0 && count <= Math.min(clean.length, offsets.length),
      "There must be one replacement offset per contaminated measurement.");
    A.check(Number.isFinite(position) && clean.every(Number.isFinite) && offsets.every(Number.isFinite),
      "Contamination inputs must be finite.");
    return clean.map((value, i) => i < count ? position + offsets[i] : value);
  }
  function compute(parameters, prepared) {
    validate(parameters);
    const data = prepared || prepare(parameters.seed);
    A.check(data.seed === parameters.seed && data.clean.length >= parameters.n && data.offsets.length >= parameters.contaminants,
      "Cached measurements do not match these settings.");
    const clean = data.clean.slice(0, parameters.n);
    const observed = contaminate(clean, data.offsets, parameters.contaminants, parameters.position);
    const influence = Array.from({ length: 121 }, (_, i) => {
      const position = -12 + i / 5;
      return { position, ...summaries(contaminate(clean, data.offsets, parameters.contaminants, position), parameters.trim) };
    });
    return { data, parameters: { ...parameters }, clean, observed, jitter: data.jitter.slice(0, parameters.n),
      summaries: summaries(observed, parameters.trim), cleanSummaries: summaries(clean, parameters.trim), influence };
  }
  return { trimmedMean, summaries, validate, prepare, contaminate, compute };
});
