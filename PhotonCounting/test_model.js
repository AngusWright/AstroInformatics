"use strict";
const assert = require("node:assert/strict");
const A = require("../TeachingTools/statistics.js"), M = require("./model.js");
const close = (actual, expected, tolerance = 1e-12) => assert.ok(Math.abs(actual - expected) <= tolerance,
  `${actual} != ${expected} (tolerance ${tolerance})`);
const originalRandom = Math.random;
Math.random = () => { throw new Error("Use the seeded RNG, not Math.random"); };
try {
  // Independent R qpois references, from tiny means through the largest supported mean.
  const probabilities = [1e-8, 0.001, 0.1, 0.5, 0.9, 0.999, 0.99999999];
  for (const [lambda, expected] of [[0, [0, 0, 0, 0, 0, 0, 0]], [0.4, [0, 0, 0, 0, 1, 3, 8]],
    [10, [0, 2, 6, 10, 14, 21, 32]], [160, [94, 122, 144, 160, 176, 200, 236]]]) {
    assert.deepEqual(probabilities.map(p => M.poissonQuantile(p, lambda)), expected);
  }
  close(M.poissonPmf(0, 0.4), 0.6703200460356393);
  close(M.poissonPmf(1, 0.4), 0.26812801841425576);
  close(M.poissonPmf(10, 10), 0.1251100357211333);
  assert.equal(M.poissonPmf(-1, 2), 0);
  assert.equal(M.poissonPmf(0, 0), 1); assert.equal(M.poissonPmf(1, 0), 0);
  assert.equal(M.poissonQuantile(1, 0), 0); assert.equal(M.poissonQuantile(1, 2), Infinity);
  assert.equal(M.poissonQuantile(0, 2), 0);
  // R pnorm(-.5,.4,sqrt(.4)); continuity-corrected Normal bin mass, not a density.
  close(A.normalCdf(-0.5, 0.4, Math.sqrt(0.4)), 0.0773644617426892);
  close(M.gaussianMass(0, 0.4), 0.562816469418554 - 0.0773644617426892);
  assert.ok(M.gaussianMass(-1, 0.4) > 0);
  assert.ok(Math.abs(M.gaussianMass(0, 0.4) - M.poissonPmf(0, 0.4)) > 0.15);
  assert.equal(M.gaussianMass(0, 0), null);

  const base = M.prepare(42), experiment = M.compute({ rate: 2.5, exposure: 4 }, base);
  // R cumsum(qexp(1-u)) for seeded uniforms .12176399340387434, .44645183894317597, .6149655369808897.
  [2.10567058806859, 2.91209433606098, 3.3982833862362].forEach((value, i) => close(base.unitArrivals[i], value));
  const zero = M.atFrame(experiment, 0), short = M.atFrame(experiment, 30), long = M.atFrame(experiment, 60);
  assert.ok(zero.counts.every(count => count === 0)); assert.deepEqual(zero.arrivals, []);
  assert.equal(zero.lambda, 0); assert.equal(zero.relativeSd, null); assert.equal(zero.gaussian, null);
  assert.equal(long.lambda, 10); close(long.sd, Math.sqrt(10)); close(long.relativeSd, 1 / Math.sqrt(10));
  close(long.mean, 10, 0.25); close(long.variance, 10, 0.7);
  assert.ok(long.arrivals.every((time, i) => time > 0 && time <= 4 && (!i || time > long.arrivals[i - 1])));
  assert.deepEqual(short.arrivals, long.arrivals.filter(time => time <= 2));
  for (let frame = 1; frame <= M.STEPS; frame++) {
    assert.ok(experiment.frames[frame].counts.every((count, i) => Number.isInteger(count) &&
      count >= experiment.frames[frame - 1].counts[i]), "Fixed-quantile coupled histogram counts must be monotone");
  }
  assert.equal(M.atFrame(experiment, 30), short, "Animation must reveal a cached frame, not sample again");
  assert.deepEqual(experiment, M.compute({ rate: 2.5, exposure: 4 }), "Seeded experiments must be repeatable");
  assert.notDeepEqual(base.uniforms, M.prepare(43).uniforms);
  const extended = M.compute({ rate: 2.5, exposure: 8 }, base);
  assert.deepEqual(extended.frames[30].arrivals, long.arrivals);
  assert.deepEqual(extended.frames[30].counts, long.counts, "Exposure changes must preserve the same underlying quantiles");
  const dark = M.compute({ rate: 0, exposure: 20 }, base);
  assert.ok(dark.frames.every(frame => frame.lambda === 0 && frame.variance === 0 &&
    frame.arrivals.length === 0 && frame.counts.every(count => count === 0) && frame.negativeMass === null));
  const low = M.compute({ rate: 0.4, exposure: 1 }, base).frames[60];
  close(low.mean, 0.4, 0.04); close(low.variance, 0.4, 0.06);
  close(low.negativeMass, 0.0773644617426892);
  close(low.empirical.reduce((sum, mass) => sum + mass, 0), 1);
  const hidden = M.compute({ rate: 2.5, exposure: 4, gaussian: "hide" }, base);
  assert.deepEqual(hidden.frames, experiment.frames, "Display choices must not change observations");
  for (const seed of [0, 2147483647]) {
    const maximum = M.compute({ seed, rate: 8, exposure: 20 }).frames[60];
    assert.equal(maximum.lambda, 160);
    close(maximum.mean, 160, 1.3); close(maximum.variance, 160, 15);
    assert.ok(maximum.counts.every(count => Number.isInteger(count) && count >= 0));
    assert.ok(maximum.exact.every(Number.isFinite) && maximum.gaussian.every(Number.isFinite));
  }
  for (const bad of [{ rate: -0.1 }, { rate: 8.1 }, { rate: NaN }, { exposure: 0 }, { exposure: 21 },
    { exposure: Infinity }, { seed: -1 }, { seed: 1.5 }, { seed: 2147483648 }, { gaussian: "approximate-sampler" }]) {
    assert.throws(() => M.compute(bad), A.ValidationError);
  }
  assert.throws(() => M.poissonQuantile(-0.1, 2), A.ValidationError);
  assert.throws(() => M.poissonQuantile(0.5, 161), A.ValidationError);
  assert.throws(() => M.poissonPmf(0.5, 2), A.ValidationError);
  assert.throws(() => M.atFrame(experiment, 61), A.ValidationError);
  console.log("PhotonCounting: R-calibrated quantiles, exact count moments, low-count Gaussian failure, cached playback and zero-rate limits passed.");
} finally { Math.random = originalRandom; }
