"use strict";
const assert = require("node:assert/strict");
const { test } = require("node:test");
const A = require("../TeachingTools/statistics.js");
const M = require("./model.js");
const base = { seed: 42, n: 200, width: 0.4, origin: 0.25, bandwidth: 0.25 };
const close = (actual, expected, tolerance = 1e-12) =>
  assert.ok(Math.abs(actual - expected) <= tolerance, `${actual} != ${expected} (tolerance ${tolerance})`);
function integrate(f, lo, hi, n = 4000) {
  const step = (hi - lo) / n;
  let sum = f(lo) + f(hi);
  for (let i = 1; i < n; i++) sum += (i % 2 ? 4 : 2) * f(lo + i * step);
  return sum * step / 3;
}

test("the known population is a normalized 60/40 two-Gaussian mixture", () => {
  close(integrate(M.populationDensity, -8, 8), 1);
  close(M.populationDensity(-0.8),
    (0.6 + 0.4 * Math.exp(-0.5 * (1.7 / 0.35) ** 2)) / (0.35 * Math.sqrt(2 * Math.PI)));
});

test("direct Gaussian KDE values and finite-window mass match R-calibrated answers", () => {
  const values = [-1, 0, 2], expected = [0.30195550199366905, 0.30204471809359046, 0.07207717550142734];
  [-1, 0, 1].forEach((x, i) => close(M.kdeAt(values, 0.5, x), expected[i]));
  close(M.kdeMass(values, 0.5, -1, 1), 0.49240606527446668);
  close(M.kdeAt([0], 0.5, 0), 2 / Math.sqrt(2 * Math.PI));
  close(integrate(x => M.kdeAt(values, 0.5, x), -10, 10), 1);
  close(M.kdeMass(values, 0.5, -20, 20), 1);
  close(M.kdeMass([0], 1, -1, 1), 0.6826894921370859);
  close(M.kdeMass([0], 1, 0, 0), 0);
  close(M.kdeOutsideMass([0], 1, -1, 1), 0.3173105078629141);
  assert.ok(M.kdeOutsideMass([0], 1, -10, 10) > 1e-23, "Tiny omitted tails must survive cancellation near one");
});

test("histogram heights are count/(n*width) and include observations on bin edges", () => {
  const data = { seed: 42, values: Array.from({ length: 50 }, (_, i) => (i % 5 - 2) / 2) };
  const result = M.compute({ ...base, n: 50, width: 0.5, origin: 0 }, data);
  assert.deepEqual(result.breaks, [-1, -0.5, 0, 0.5, 1, 1.5]);
  assert.deepEqual(result.histogram.counts, [10, 10, 10, 10, 10]);
  result.histogram.density.forEach(height => close(height, 0.4));
  close(result.histogramArea, 1);
  assert.equal(result.histogram.excluded, 0);
  const shift = M.compute({ ...base, n: 50, width: 0.5, origin: 1 }, data);
  assert.deepEqual(shift.histogram, result.histogram, "One whole bin shift is the same partition");
  const decimalEdges = { seed: 42, values: Array(50).fill(0.3) };
  const decimal = M.compute({ ...base, n: 50, width: 0.1, origin: 0 }, decimalEdges);
  close(decimal.histogramArea, 1);
  assert.equal(decimal.histogram.excluded, 0);
});

test("all smoothing controls preserve the exact observations while changing their display", () => {
  const data = M.prepare(42), first = M.compute(base, data);
  const moved = M.compute({ ...base, origin: 0.75 }, data);
  const wider = M.compute({ ...base, width: 0.8 }, data);
  const smoother = M.compute({ ...base, bandwidth: 1 }, data);
  for (const result of [moved, wider, smoother]) {
    assert.strictEqual(result.data, first.data);
    assert.deepEqual(result.values, first.values);
    close(result.histogramArea, 1);
  }
  assert.notDeepEqual(first.histogram.counts, moved.histogram.counts);
  assert.notDeepEqual(first.histogram.counts, wider.histogram.counts);
  assert.deepEqual(first.density, moved.density);
  assert.deepEqual(first.density, wider.density);
  assert.notDeepEqual(first.density, smoother.density);
  assert.ok(smoother.outsideMass > 0.01, "Broad Gaussian tails must not be silently renormalised away");
  close(smoother.outsideMass + M.kdeMass(smoother.values, 1, ...smoother.limits), 1);
  const larger = M.compute({ ...base, n: 1000 }, data);
  assert.deepEqual(first.values, larger.values.slice(0, base.n));
  assert.deepEqual(M.prepare(42), data);
});

test("bandwidth reveals sampling bumps or merges real populations, not a component-count test", () => {
  const data = M.prepare(42);
  const sharp = M.compute({ ...base, n: 1000, bandwidth: 0.15 }, data);
  const smooth = M.compute({ ...base, n: 1000, bandwidth: 1 }, data);
  assert.ok(sharp.peaks.length >= 2);
  assert.equal(smooth.peaks.length, 1);
  assert.deepEqual(M.countPeaks([0, 0.01, 0, 1, 0]), [3]);
  assert.deepEqual(M.countPeaks([0, 1, 1, 2, 0]), [3], "A shoulder is not an extra peak");
  assert.deepEqual(M.countPeaks([0, 2, 2, 0]), [1], "A flat-topped maximum counts once");
});

test("invalid smoothing, sample settings, KDE arguments and caches are rejected", () => {
  for (const patch of [{ seed: -1 }, { seed: 0.5 }, { n: 49 }, { n: 1001 }, { n: 50.5 },
    { width: 0.09 }, { width: 1.51 }, { origin: -0.1 }, { origin: 1.1 }, { bandwidth: 0 }, { bandwidth: 1.1 }]) {
    assert.throws(() => M.compute({ ...base, ...patch }), A.ValidationError);
  }
  for (const [values, bandwidth] of [[[], 1], [[NaN], 1], [[0], 0], [[0], Infinity]]) {
    assert.throws(() => M.kdeAt(values, bandwidth, 0), A.ValidationError);
  }
  assert.throws(() => M.kdeMass([0], 1, 2, 1), A.ValidationError);
  assert.throws(() => M.countPeaks([0, -1, 0]), A.ValidationError);
  assert.throws(() => M.compute(base, M.prepare(0)), A.ValidationError);
});
