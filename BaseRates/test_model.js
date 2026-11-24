"use strict";
const assert = require("node:assert/strict");
const A = require("../TeachingTools/statistics.js");
const M = require("./model.js");
const close = (actual, expected) => assert.ok(Math.abs(actual - expected) < 1e-13, `${actual} != ${expected}`);

// Retained model.R references and independent Bayes-rule natural frequencies.
close(M.positivePredictiveValue(0.01, 0.95, 0.01), 0.48969072164948452);
close(M.positivePredictiveValue(0.0001, 0.95, 0.01), 0.0094115316029324356);
close(M.positivePredictiveValue(0.001, 0.95, 0.01), 950 / (950 + 9990));
assert.ok(M.positivePredictiveValue(0.001, 0.95, 0.001) > M.positivePredictiveValue(0.001, 0.95, 0.01));
assert.equal(M.positivePredictiveValue(0.001, 0.95, 0), 1);
assert.equal(M.positivePredictiveValue(0, 0.95, 0), null);
assert.equal(M.positivePredictiveValue(1, 0, 0.1), null);
const expected = M.expectedCounts(10000, 0.0001, 0.95, 0.01);
close(expected.truePositive, 0.95);
close(expected.falsePositive, 99.99);
close(Object.values(expected).reduce((a, b) => a + b, 0), 10000);

const p = { seed: 42, log_prevalence: -2, sensitivity: 0.95, false_positive: 0.01, reveal: "show" };
const survey = M.compute(p), hidden = M.compute({ ...p, reveal: "hide" });
assert.deepEqual(survey.truth, hidden.truth);
assert.deepEqual(survey.flagged, hidden.flagged);
assert.deepEqual(survey, M.compute(p));
assert.notDeepEqual(survey.truth, M.compute({ ...p, seed: 43 }).truth);
const hand = { truth: [true, false, true, false], flagged: [true, true, false, false],
  prevalence: 0.5, sensitivity: 0.5, falsePositive: 0.5 };
assert.deepEqual(M.reveal(hand).actual, { truePositive: 1, falsePositive: 1, falseNegative: 1, trueNegative: 1 });
assert.equal(M.reveal(hand).genuineFraction, 0.5);
assert.equal(M.reveal(hand, 0).genuineFraction, null);
assert.equal(M.reveal(hand, 50).indices.length, 2);
const empty = M.compute({ ...p, log_prevalence: -6, false_positive: 0 });
assert.equal(M.reveal(empty).indices.length, 0);
assert.equal(M.reveal(empty).genuineFraction, null);
assert.equal(empty.probability, 1, "An empty observed survey must not turn theoretical PPV into zero");
const rare = M.compute({ ...p, log_prevalence: -6 });
assert.ok(rare.probability > 0 && rare.probability < 0.001);
assert.ok(M.reveal(rare).expected.truePositive > 0 && M.reveal(rare).expected.truePositive < 1);
assert.equal(M.reveal(survey, 41).scanned, 4100);
assert.deepEqual(M.reveal(survey, 41).indices, M.reveal(survey).indices.filter(i => i < 4100));
for (const bad of [{ seed: -1 }, { log_prevalence: -7 }, { sensitivity: 1.1 }, { false_positive: -0.01 }, { reveal: "invalid" }])
  assert.throws(() => M.compute({ ...p, ...bad }), A.ValidationError);
assert.throws(() => M.positivePredictiveValue(NaN, 0.9, 0.1), A.ValidationError);
assert.throws(() => M.reveal(survey, -1), A.ValidationError);
console.log("BaseRates: R/Bayes references, rare fractional expectations, empty flags and fixed-identity replay passed.");
