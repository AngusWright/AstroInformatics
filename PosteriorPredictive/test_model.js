"use strict";
const assert = require("node:assert/strict");
const A = require("../TeachingTools/statistics.js");
const M = require("./model.js");
const close = (actual, expected, tolerance = 1e-12) => assert.ok(Math.abs(actual - expected) < tolerance, `${actual} != ${expected}`);

// Exact retained test_numerics.R example; variances are not standard deviations.
const oracle = M.predictiveNormal([1, 2, 3], 1, 2);
close(oracle.mean, 1.875);
assert.deepEqual(Object.values(oracle.variance), [1.5625, 2.5625, 6.5625]);
const reference = M.predictiveNormal([1, 2, 3], 1.5, 0.8);
close(reference.mean, 1.9257927846963665);
close(reference.variance.mean, 0.92759019129541664);
close(reference.variance.star, 3.17759019129541675);
close(reference.variance.measurement, 3.81759019129541688);
// R qnorm(.975) endpoints, not the rounded 1.96 approximation.
close(reference.intervals.measurement[0], -1.9037131252281956);
close(reference.intervals.measurement[1], 5.7552986946209286);
const large = M.predictiveNormal(Array(100000).fill(2), 1, 2);
assert.ok(large.variance.mean < 0.0001);
close(large.variance.star - large.variance.mean, 1);
close(large.variance.measurement - large.variance.star, 4);
const zeroScatter = M.predictiveNormal([2], 0, 0.1);
close(zeroScatter.mean, 1.9992003198720509);
assert.equal(zeroScatter.variance.mean, zeroScatter.variance.star);
assert.deepEqual(zeroScatter.intervals.mean, zeroScatter.intervals.star);
close(zeroScatter.variance.measurement - zeroScatter.variance.star, 0.01);

const p = { seed: 42, n: 20, intrinsic: 1.5, measurement: 0.8, show_data: true };
const first = M.compute(p), many = M.compute({ ...p, n: 300 }), hidden = M.compute({ ...p, show_data: false });
assert.deepEqual(first.x, many.x.slice(0, 20));
assert.deepEqual(first.x, hidden.x);
assert.deepEqual(first.observations, many.observations);
assert.deepEqual(first, M.compute(p));
assert.notDeepEqual(first.x, M.compute({ ...p, seed: 43 }).x);
assert.ok(many.posterior.variance.mean < first.posterior.variance.mean);
assert.ok(many.posterior.variance.measurement > 2.89, "More stars cannot remove intrinsic and measurement variance");
assert.equal(M.compute({ ...p, intrinsic: 0 }).stars.every(x => x === 2), true);
M.compute({ ...p, seed: 2147483647, n: 1, intrinsic: 0, measurement: 0.1 });
for (const bad of [{ seed: NaN }, { n: 0 }, { n: 1.5 }, { intrinsic: -1 }, { measurement: 0 }, { show_data: "yes" }])
  assert.throws(() => M.compute({ ...p, ...bad }), A.ValidationError);
assert.throws(() => M.predictiveNormal([], 1, 2), A.ValidationError);
assert.throws(() => M.predictiveNormal([1, Infinity], 1, 2), A.ValidationError);
assert.throws(() => M.predictiveNormal([1], 1, 2, 0), A.ValidationError);
console.log("PosteriorPredictive: retained R variance oracle, qnorm intervals, irreducible floors and nested observations passed.");
