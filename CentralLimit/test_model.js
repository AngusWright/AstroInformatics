"use strict";
const assert = require("node:assert/strict");
const { test } = require("node:test");
const A = require("../TeachingTools/statistics.js");
const M = require("./model.js");
const base = { seed: 42, population: "exponential", n: 100, repetitions: 3000 };
const close = (actual, expected, tolerance = 1e-12) =>
  assert.ok(Math.abs(actual - expected) <= tolerance, `${actual} != ${expected} (tolerance ${tolerance})`);

test("analytic population moments and densities are not fitted from a finite run", () => {
  close(M.populations.exponential.mean, 1);
  close(M.populations.exponential.variance, 1);
  close(M.populations.mixture.mean, 0.05 * 5);
  close(M.populations.mixture.variance, 1 + 0.05 * 0.95 * 25);
  assert.equal(M.populations.cauchy.mean, null);
  assert.equal(M.populations.cauchy.variance, null);
  close(M.populationDensity("exponential", 1), Math.exp(-1));
  assert.equal(M.populationDensity("exponential", -0.01), 0);
  close(M.populationDensity("mixture", 0), (0.95 + 0.05 * Math.exp(-12.5)) / Math.sqrt(2 * Math.PI));
  close(M.populationDensity("cauchy", 1), 1 / (2 * Math.PI));
});

test("frame arithmetic uses a fixed independent sequence for each experiment", () => {
  const data = { seed: 42, population: "exponential", observations: [], sums: new Float64Array(300000) };
  for (let i = 0; i < 200; i++) {
    for (let n = 1; n <= 100; n++) data.sums[i * 100 + n - 1] = n * (n + 1) / 2 + (i % 2) * n;
  }
  const result = M.compute({ ...base, repetitions: 200 }, data);
  const start = M.atFrame(result, 0), middle = M.atFrame(result, 50), finish = M.atFrame(result);
  assert.equal(start.n, 1); assert.equal(middle.n, 51); assert.equal(finish.n, 100);
  assert.deepEqual(start.means.slice(0, 4), [1, 2, 1, 2]);
  assert.deepEqual(middle.means.slice(0, 2), [26, 27]);
  assert.deepEqual(finish.means.slice(0, 2), [50.5, 51.5]);
  assert.deepEqual(finish.statistic.slice(0, 2), [495, 505]);
  close(finish.meanSd, 0.1);
  close(finish.empiricalIqr, 1);
  assert.deepEqual(M.atFrame(result), finish, "Playback cannot draw new observations");
});

test("finite-variance sample means converge with the correct mixture variance", () => {
  for (const population of ["exponential", "mixture"]) {
    const prepared = M.prepare(42, population), result = M.compute({ ...base, population }, prepared);
    const finish = M.atFrame(result), summary = A.summarize(finish.statistic);
    assert.ok(Math.abs(summary.mean) < 0.08, `${population}: standardised mean ${summary.mean}`);
    assert.ok(Math.abs(summary.sd ** 2 - 1) < 0.1, `${population}: standardised variance ${summary.sd ** 2}`);
    close(finish.meanSd, population === "mixture" ? Math.sqrt(2.1875) / 10 : 0.1);
    const fewer = M.compute({ ...base, population, repetitions: 200 }, prepared);
    assert.strictEqual(fewer.data, result.data);
    assert.deepEqual(M.atFrame(fewer).means, finish.means.slice(0, 200));
    const one = M.compute({ ...base, population, n: 1 }, prepared);
    assert.strictEqual(one.data, result.data);
    assert.deepEqual(M.atFrame(one).means, M.atFrame(result, 0).means);
  }
});

test("Cauchy averages keep their exact population law, not a fictional finite SD", () => {
  const result = M.compute({ ...base, population: "cauchy" });
  const start = M.atFrame(result, 0), finish = M.atFrame(result);
  assert.equal(finish.finite, false);
  assert.equal(finish.meanSd, null);
  assert.strictEqual(finish.means, finish.statistic);
  for (const run of [start, finish]) {
    close(A.quantile(run.means, 0.25), -1, 0.16);
    close(A.quantile(run.means, 0.5), 0, 0.1);
    close(A.quantile(run.means, 0.75), 1, 0.16);
    close(run.empiricalIqr, 2, 0.3);
  }
  assert.ok(finish.empiricalIqr / start.empiricalIqr > 0.7, "Cauchy widths must not shrink as 1/sqrt(n)");
});

test("out-of-window Cauchy tails remain in the histogram denominator", () => {
  const hist = A.histogram([-100, -1, 0, 1, 100], [-10, 0, 10]);
  assert.equal(hist.excluded, 2);
  close(hist.density.reduce((sum, height, i) => sum + height * hist.widths[i], 0), 3 / 5);
});

test("invalid settings and mismatched caches are expected validation errors", () => {
  for (const patch of [{ seed: -1 }, { seed: 0.5 }, { n: 0 }, { n: 101 }, { n: 1.5 },
    { repetitions: 199 }, { repetitions: 3001 }, { population: "normal" }, { population: "toString" }]) {
    assert.throws(() => M.compute({ ...base, ...patch }), A.ValidationError);
  }
  const prepared = M.prepare(0, "exponential");
  assert.throws(() => M.compute(base, prepared), A.ValidationError);
  const result = M.compute({ ...base, seed: 0 }, prepared);
  for (const frame of [-1, 101, 0.5, NaN]) assert.throws(() => M.atFrame(result, frame), A.ValidationError);
});
