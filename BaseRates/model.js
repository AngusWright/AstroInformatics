(function (root) {
  "use strict";
  const A = typeof module === "object" && module.exports ? require("../TeachingTools/statistics.js") : root.Astro;
  function probability(value, label) {
    A.check(Number.isFinite(value) && value >= 0 && value <= 1, `${label} must be between 0 and 1.`);
  }
  function positivePredictiveValue(prevalence, sensitivity, falsePositive) {
    probability(prevalence, "Prevalence"); probability(sensitivity, "Sensitivity"); probability(falsePositive, "False-positive rate");
    const real = prevalence * sensitivity, falseFlag = (1 - prevalence) * falsePositive;
    return real + falseFlag === 0 ? null : real / (real + falseFlag);
  }
  function expectedCounts(n, prevalence, sensitivity, falsePositive) {
    A.check(Number.isInteger(n) && n >= 0, "Survey size must be a non-negative integer.");
    positivePredictiveValue(prevalence, sensitivity, falsePositive);
    return { truePositive: n * prevalence * sensitivity, falsePositive: n * (1 - prevalence) * falsePositive,
      falseNegative: n * prevalence * (1 - sensitivity), trueNegative: n * (1 - prevalence) * (1 - falsePositive) };
  }
  function compute(p) {
    A.check(Number.isInteger(p.seed) && p.seed >= 0 && p.seed <= 2147483647, "Seed must be an integer from 0 to 2147483647.");
    A.check(Number.isFinite(p.log_prevalence) && p.log_prevalence >= -6 && p.log_prevalence <= -1, "log10 prevalence must be between -6 and -1.");
    A.check(Number.isFinite(p.sensitivity) && p.sensitivity >= 0.5 && p.sensitivity <= 1, "Sensitivity must be between 0.5 and 1.");
    A.check(Number.isFinite(p.false_positive) && p.false_positive >= 0 && p.false_positive <= 0.1, "False-positive rate must be between 0 and 0.1.");
    A.check(["show", "hide"].includes(p.reveal), "Choose whether to reveal the truth.");
    const prevalence = 10 ** p.log_prevalence, rng = A.createRng(p.seed);
    const truth = Array.from({ length: 10000 }, () => rng.uniform() < prevalence);
    const flagged = truth.map(real => rng.uniform() < (real ? p.sensitivity : p.false_positive));
    return { truth, flagged, prevalence, sensitivity: p.sensitivity, falsePositive: p.false_positive,
      probability: positivePredictiveValue(prevalence, p.sensitivity, p.false_positive) };
  }
  function reveal(result, frame = 100) {
    A.check(Number.isFinite(frame) && frame >= 0 && frame <= 100, "Survey progress must be between 0 and 100.");
    const scanned = Math.floor(result.truth.length * frame / 100), indices = [];
    const actual = { truePositive: 0, falsePositive: 0, falseNegative: 0, trueNegative: 0 };
    for (let i = 0; i < scanned; i++) {
      const real = result.truth[i], flag = result.flagged[i];
      actual[flag ? (real ? "truePositive" : "falsePositive") : (real ? "falseNegative" : "trueNegative")]++;
      if (flag) indices.push(i);
    }
    return { scanned, indices, actual, expected: expectedCounts(scanned, result.prevalence, result.sensitivity, result.falsePositive),
      genuineFraction: indices.length ? actual.truePositive / indices.length : null };
  }
  const api = { positivePredictiveValue, expectedCounts, compute, reveal };
  if (typeof module === "object" && module.exports) module.exports = api;
  else root.BaseRatesModel = api;
})(globalThis);
