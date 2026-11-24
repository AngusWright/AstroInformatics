"use strict";
const assert = require("node:assert/strict");
const { test } = require("node:test");
const A = require("../TeachingTools/statistics.js");
const M = require("./model.js");
const base = { seed: 42, family: "power", slope: 2.35, n: 500 };
const close = (actual, expected, tolerance = 1e-12) =>
  assert.ok(Math.abs(actual - expected) <= tolerance, `${actual} != ${expected} (tolerance ${tolerance})`);
function integrate(f, lo, hi, n = 4000) {
  const step = (hi - lo) / n;
  let sum = f(lo) + f(hi);
  for (let i = 1; i < n; i++) sum += (i % 2 ? 4 : 2) * f(lo + i * step);
  return sum * step / 3;
}

test("bounded power law agrees with independent R-calibrated values", () => {
  const target = M.distribution("power", 2.35);
  const expectedQuantiles = [0.10386943771591439, 0.16700653978735835, 0.90986487886954637];
  [0.05, 0.5, 0.95].forEach((p, i) => close(target.quantile(p), expectedQuantiles[i]));
  [0.1, 1, 20].forEach((x, i) => close(target.density(x), [13.510575072192447, 0.060349522052721795, 0.000052875360962235123][i]));
  [0.2, 1, 10].forEach((x, i) => close(target.cdf(x), [0.6081839921963873, 0.9560799890122392, 0.9987865133998033][i]));
});

test("alpha=1 is the exact normalized logarithmic limit, including endpoints", () => {
  const logRange = Math.log(200);
  for (const slope of [0.5, 1 - 1e-10, 1, 1 + 1e-10, 2.35, 3.5]) {
    const target = M.distribution("power", slope);
    assert.equal(target.quantile(0), 0.1);
    assert.equal(target.quantile(1), 20);
    assert.equal(target.cdf(-1), 0); assert.equal(target.cdf(0.1), 0);
    assert.equal(target.cdf(20), 1); assert.equal(target.cdf(Infinity), 1);
    assert.equal(target.density(0), 0); assert.equal(target.density(20.1), 0);
    for (const p of [1e-8, 0.05, 0.5, 0.99, 1 - 1e-8]) close(target.cdf(target.quantile(p)), p, 2e-13);
    close(integrate(z => M.logDensity(target, z), -1, Math.log10(20)), 1, 2e-10);
    if (Math.abs(slope - 1) < 1e-8) {
      close(target.quantile(0.5), Math.sqrt(2), 1e-9);
      close(target.density(1), 1 / logRange, 1e-10);
    }
  }
  close(M.logDensity(M.distribution("power", 1), 0), 1 / Math.log10(200));
});

test("exponential and log-normal retain normalized targets and valid endpoints", () => {
  const exponential = M.distribution("exponential"), lognormal = M.distribution("lognormal");
  close(exponential.quantile(0.5), Math.log(2));
  close(exponential.cdf(Math.log(2)), 0.5);
  close(lognormal.quantile(0.5), 1);
  close(lognormal.cdf(1), 0.5);
  close(lognormal.density(1), 1 / (0.7 * Math.sqrt(2 * Math.PI)));
  for (const target of [exponential, lognormal]) {
    assert.equal(target.quantile(0), 0); assert.equal(target.quantile(1), Infinity);
    assert.equal(target.cdf(0), 0); assert.equal(target.cdf(Infinity), 1);
    assert.equal(target.density(Infinity), 0);
    close(integrate(z => M.logDensity(target, z), -12, 3), 1, 2e-10);
    for (const p of [0.00001, 0.05, 0.5, 0.9, 0.99999]) close(target.cdf(target.quantile(p)), p, 2e-12);
  }
});

test("the two-sided KS discrepancy includes both sides of each empirical jump", () => {
  close(M.ksStatistic([0.1, 0.4, 0.6, 0.9]), 0.15);
  close(M.ksStatistic([0.8, 0.9]), 0.8);
  close(M.ksStatistic([0, 0, 1]), 2 / 3);
  assert.equal(M.ksStatistic([]), null);
});

test("family, sample size and animation reuse one cached uniform sequence", () => {
  const prepared = M.prepare(42), power = M.compute(base, prepared);
  const exponential = M.compute({ ...base, family: "exponential", n: 1500 }, prepared);
  assert.strictEqual(power.data, exponential.data);
  assert.deepEqual(power.uniforms, exponential.uniforms.slice(0, 500));
  assert.notDeepEqual(power.values, exponential.values.slice(0, 500));
  assert.equal(M.atFrame(power, 0).discrepancy, null);
  assert.equal(M.atFrame(power, 0).count, 0);
  assert.equal(M.atFrame(power, 1).count, 5);
  assert.deepEqual(M.atFrame(power, 50).values, power.values.slice(0, 250));
  assert.deepEqual(M.atFrame(power).values, power.values);
  assert.deepEqual(M.prepare(42).uniforms, prepared.uniforms);
  assert.notDeepEqual(M.prepare(43).uniforms, prepared.uniforms);
  for (const result of [power, exponential]) {
    assert.ok(result.values.every(x => Math.log10(x) >= result.limits[0] && Math.log10(x) <= result.limits[1]));
  }
});

test("invalid probabilities, parameters, frames and caches are rejected", () => {
  for (const family of ["power", "exponential", "lognormal"]) {
    for (const p of [-0.01, 1.01, NaN, Infinity]) assert.throws(() => M.distribution(family).quantile(p), A.ValidationError);
  }
  for (const patch of [{ seed: -1 }, { seed: 0.5 }, { n: 49 }, { n: 1501 }, { n: 50.5 },
    { family: "unknown" }, { slope: 0.49 }, { slope: 3.51 }, { slope: NaN }]) {
    assert.throws(() => M.compute({ ...base, ...patch }), A.ValidationError);
  }
  assert.throws(() => M.compute(base, M.prepare(0)), A.ValidationError);
  assert.throws(() => M.ksStatistic([1.1]), A.ValidationError);
  assert.throws(() => M.atFrame(M.compute(base), -1), A.ValidationError);
});
