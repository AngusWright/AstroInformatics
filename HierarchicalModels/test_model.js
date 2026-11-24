"use strict";
const assert = require("node:assert/strict");
const { test } = require("node:test");
const M = require("./model.js");
const { ValidationError } = require("../TeachingTools/statistics.js");
const close = (actual, expected, tolerance = 2e-12) => assert.ok(Math.abs(actual - expected) < tolerance,
  `${actual} != ${expected}`);

test("normal-normal inference agrees with independent R posterior values", () => {
  const p = M.normalHierarchy([0, 2, 5], [0.25, 1, 2], 1);
  [0.069979006298110588, 1.59482155353393984, 1.951714485654303832]
    .forEach((x, i) => close(p.mean[i], x));
  [0.060881735479356193, 0.648705388383485015, 1.180685794261721711]
    .forEach((x, i) => close(p.variance[i], x));
  close(p.populationMean, 1.1896431070678797);
  close(p.populationSd, 0.77124675268939691);
  [16 / 17, 0.5, 0.2].forEach((x, i) => close(p.weight[i], x));
  assert.ok(p.variance[2] > 0.8, "Uncertainty in the shared mean must be propagated, not plugged in.");
});

test("marginal intervals match a separate multivariate Gaussian covariance calculation", () => {
  // R oracle: C0 = diag(0.7^2,3) + matrix(25,3,3);
  // C = solve(solve(C0) + diag(1/c(.3,.8,2)^2)); mean = C %*% (c(-1,1,2)/c(.3,.8,2)^2).
  const p = M.normalHierarchy([-1, 1, 2], [0.3, 0.8, 2], 0.7);
  [-0.86610294610661565, 0.3559743009363221, 0.096117706094160815]
    .forEach((value, i) => close(p.mean[i], value));
  [0.084418906178909389, 0.3892205957720189, 0.71288324993832797]
    .forEach((value, i) => close(p.variance[i], value));
});

test("zero scatter is uncertain complete pooling, not zero-width intervals", () => {
  const p = M.normalHierarchy([0, 2, 5], [0.25, 1, 2], 0);
  p.mean.forEach(x => close(x, 0.18796992481203009));
  p.variance.forEach(x => close(x, 0.057836899942163102));
  assert.deepEqual(p.weight, [0, 0, 0]);
  const shiftedPrior = M.normalHierarchy([2], [1], 0, 10, 2);
  close(shiftedPrior.populationMean, 3.6);
  close(shiftedPrior.variance[0], 0.8);
  const partial = M.compute({ scatter: 0 }), complete = M.compute({ scatter: 0, pooling: "complete" });
  assert.deepEqual(partial.posterior, complete.posterior);
  const independent = M.normalHierarchy([-1, 1, 2], [0.3, 0.8, 2], 1e6);
  independent.mean.forEach((value, i) => close(value, [-1, 1, 2][i], 1e-8));
  independent.variance.forEach((value, i) => close(value, [0.09, 0.64, 4][i], 1e-8));
});

test("pooling comparisons hold data fixed and interpolation reuses the inference", () => {
  const partial = M.compute(), none = M.compute({ pooling: "none" }), complete = M.compute({ pooling: "complete" });
  assert.deepEqual(partial, M.compute());
  for (const other of [none, complete]) {
    assert.deepEqual(partial.observed, other.observed);
    assert.deepEqual(partial.truth, other.truth);
  }
  assert.deepEqual(none.final.estimate, none.observed);
  none.final.sd.forEach((x, i) => close(x, none.error[i]));
  assert.equal(none.posterior.populationMean, null);
  assert.ok(complete.posterior.mean.every(x => x === complete.posterior.populationMean));
  assert.ok(complete.posterior.variance.every(x => x > 0));
  assert.ok(partial.posterior.weight[0] > partial.posterior.weight[13]);
  assert.deepEqual(M.atPhase(partial, 0).estimate, partial.observed);
  assert.deepEqual(M.atPhase(partial, 0).sd, partial.error);
  const middle = M.atPhase(partial, 0.5);
  middle.estimate.forEach((x, i) => close(x, (partial.observed[i] + partial.posterior.mean[i]) / 2));
  assert.deepEqual(M.atPhase(partial), partial.final);
  assert.notDeepEqual(partial.observed, M.compute({ seed: 43 }).observed);
});

test("pure API validates domains and supports all UI dispersion/noise boundaries", () => {
  for (const scatter of [0, 3]) for (const noise of [0.3, 3]) {
    assert.ok(M.compute({ scatter, noise }).posterior.variance.every(x => Number.isFinite(x) && x > 0));
  }
  for (const p of [{ scatter: -1 }, { noise: 0 }, { pooling: "bad" }, { seed: 0.5 }]) {
    assert.throws(() => M.compute(p), ValidationError);
  }
  assert.throws(() => M.normalHierarchy([], [], 1), ValidationError);
  assert.throws(() => M.normalHierarchy([1], [0], 1), ValidationError);
  assert.throws(() => M.normalHierarchy([1], [1], 1, 0, Infinity), ValidationError);
  assert.throws(() => M.atPhase(M.compute(), NaN), ValidationError);
});
