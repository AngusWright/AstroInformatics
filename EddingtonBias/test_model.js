"use strict";
const assert = require("node:assert/strict");
const { test } = require("node:test");
const M = require("./model.js");
const { ValidationError } = require("../TeachingTools/statistics.js");
const close = (actual, expected, tolerance = 2e-12) => assert.ok(Math.abs(actual - expected) < tolerance,
  `${actual} != ${expected}`);

test("finite-boundary convolution agrees with independently evaluated R fixtures", () => {
  const xs = [7.5, 8.5, 9, 10.5, 12, 12.5];
  const expected = [0.00057858657576764609, 0.54455168729052783, 0.70001348641446304,
    0.083073993335735438, 0.0058983897032223696, 0.00046267704137797775];
  xs.forEach((x, i) => close(M.massDensity(x, 1.5, 0.3), expected[i]));
  [0.57699223421809442, 0.62790741039772291, 0.00077398456401489709]
    .forEach((value, i) => close(M.massDensity([8.5, 9, 12][i], 2.5, 0.5), value));
  assert.equal(M.massDensity(8.49, 1.5), 0);
  assert.equal(M.massDensity(12.01, 1.5), 0);
  close(M.massDensity(8.5, 1.5), 1.5 / (1 - Math.exp(-1.5 * 3.5)));
  assert.ok(M.massDensity(8.3, 1.5, 0.02) > 0, "Positive-tail CDF cancellation must not erase a resolvable density.");
});

test("both analytic densities retain unit area, including Gaussian tails", () => {
  function integral(fn, lo, hi, n = 12000) {
    const h = (hi - lo) / n;
    let sum = fn(lo) + fn(hi);
    for (let i = 1; i < n; i++) sum += (i % 2 ? 4 : 2) * fn(lo + i * h);
    return sum * h / 3;
  }
  close(integral(x => M.massDensity(x, 3), 8.5, 12), 1, 1e-10);
  close(integral(x => M.massDensity(x, 3, 0.6), 3.5, 17), 1, 1e-10);
  close(integral(x => M.massDensity(x, 0.3, 0.02), 7.5, 13), 1, 1e-10);
});

test("zero noise, cached errors and all migration counts conserve the population", () => {
  const zero = M.compute({ noise: 0 }), noisy = M.compute({ noise: 0.6 });
  assert.deepEqual(noisy, M.compute({ noise: 0.6 }));
  assert.deepEqual(zero.truth, noisy.truth);
  assert.deepEqual(zero.errors, noisy.errors);
  assert.deepEqual(zero.final.observed, zero.truth);
  assert.equal(zero.final.up, 0); assert.equal(zero.final.down, 0);
  assert.ok(noisy.truth.every(x => x >= 8.5 && x <= 12));
  for (const phase of [0, 0.25, 0.7, 1]) {
    const s = M.atPhase(noisy, phase);
    assert.equal(s.observedHigh, s.trueHigh + s.up - s.down);
    assert.equal(s.lowLow + s.up + s.down + s.highHigh, 8000);
    assert.equal(s.histogram.counts.reduce((a, b) => a + b, 0), 8000);
    assert.equal(s.histogram.excluded, 0);
    s.observed.forEach((x, i) => close(x, noisy.truth[i] + phase * 0.6 * noisy.errors[i]));
  }
  assert.deepEqual(noisy.final.observed, M.compute({ noise: 0.6, threshold: 11.5 }).final.observed);
  assert.notDeepEqual(noisy.truth, M.compute({ seed: 43 }).truth);
});

test("boundary controls remain valid and invalid pure calls fail explicitly", () => {
  for (const slope of [0.3, 3]) for (const threshold of [9.5, 11.5]) {
    assert.equal(M.compute({ noise: 0.6, slope, threshold }).final.histogram.excluded, 0);
  }
  for (const p of [{ noise: -0.1 }, { noise: 0.7 }, { slope: 0 }, { threshold: 12 }, { seed: -1 }]) {
    assert.throws(() => M.compute(p), ValidationError);
  }
  assert.throws(() => M.massDensity(9, 0), ValidationError);
  assert.throws(() => M.atPhase(M.compute(), -1), ValidationError);
});
