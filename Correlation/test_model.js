"use strict";
const assert = require("node:assert/strict");
const A = require("../TeachingTools/statistics.js"), M = require("./model.js");
const close = (actual, expected, tolerance = 1e-12) => assert.ok(Math.abs(actual - expected) <= tolerance,
  `${actual} != ${expected} (tolerance ${tolerance})`);
const originalRandom = Math.random;
Math.random = () => { throw new Error("Use the seeded RNG, not Math.random"); };
try {
  const known = M.statistics([1, 2, 3, 4], [1, 2, 1, 2]);
  close(known.correlation, 1 / Math.sqrt(5)); close(known.slope, 0.2); close(known.intercept, 1);
  close(M.statistics([1, 2, 3], [3, 2, 1]).correlation, -1);
  assert.equal(M.statistics([1, 1, 1], [1, 2, 3]).correlation, null);
  assert.equal(M.statistics([1, 2, 3], [1, 1, 1]).correlation, null);
  assert.equal(M.statistics([], []).correlation, null);
  assert.equal(M.statistics([1, 2], [3, 4]).correlation, null, "A highlighted slice needs at least three objects");
  close(M.statistics([0, 1e-140, 2e-140], [0, 2e-140, 4e-140]).correlation, 1);
  close(M.statistics([1e10 + 1, 1e10 + 2, 1e10 + 3], [-3, -2, -1]).correlation, 1);

  // Z~U(0,1), hence Cov(X,Y)=4 Var(Z)=1/3, while each variance is 1/3 + sigma^2.
  close(M.populationCorrelation(0.4), 25 / 37);
  close(M.populationCorrelation(0.4, 0.2), 1 / 13);
  close(M.populationCorrelation(2), 1 / 13);
  close(M.populationCorrelation(0), 1);
  close(M.populationCorrelation(0.4, 0), 0);
  assert.equal(M.populationCorrelation(0, 0), null);
  const base = M.prepare(42), full = M.compute({ width: 1 }, base);
  const narrow = M.compute({ width: 0.05 }, base), edge = M.compute({ centre: 0, width: 0.2 }, base);
  assert.equal(full.slice.n, M.SIZE); close(full.full.correlation, full.slice.correlation);
  close(full.full.correlation, 25 / 37, 0.05);
  close(full.intrinsic.correlation, 0, 0.08);
  assert.ok(Math.abs(narrow.slice.correlation) < 0.3 && Math.abs(narrow.slice.correlation) < full.full.correlation);
  assert.ok(narrow.slice.n > 3 && narrow.slice.n < full.slice.n);
  assert.equal(edge.lower, 0); close(edge.upper, 0.1); close(edge.theorySlice, 1 / 49);
  assert.ok(edge.selected.every((selected, i) => selected === (edge.distance[i] <= 0.1)));
  assert.deepEqual(full.x, narrow.x); assert.deepEqual(full.y, narrow.y);
  assert.equal(full.distance, narrow.distance, "Display-only slice changes must retain the population");
  assert.deepEqual(full, M.compute({ width: 1 }), "Seeded population must be repeatable");
  assert.notDeepEqual(full.distance, M.compute({ seed: 43 }).distance);
  const rescaled = M.compute({ scatter: 0.8 }, base);
  rescaled.intrinsicX.forEach((value, i) => close(value, 2 * full.intrinsicX[i]));
  const zero = M.compute({ scatter: 0 }, base);
  close(zero.full.correlation, 1); close(zero.slice.correlation, 1);
  assert.equal(zero.intrinsic.correlation, null, "Degenerate intrinsic properties have undefined Pearson correlation");
  assert.ok(zero.intrinsicX.every(value => value === 0));
  for (const seed of [0, 2147483647]) for (const scatter of [0, 2]) for (const centre of [0, 1]) {
    const result = M.compute({ seed, scatter, centre, width: 0.05 });
    assert.ok(result.x.every(Number.isFinite) && result.y.every(Number.isFinite));
    assert.ok(result.slice.n >= 3);
    if (scatter === 0) assert.equal(result.intrinsic.correlation, null);
  }
  for (const bad of [{ scatter: -0.1 }, { scatter: NaN }, { scatter: 2.1 }, { centre: -0.01 }, { centre: 1.1 },
    { width: 0 }, { width: 1.1 }, { seed: -1 }, { seed: 0.5 }, { seed: 2147483648 }]) {
    assert.throws(() => M.compute(bad), A.ValidationError);
  }
  assert.throws(() => M.statistics([1, 2, 3], [1, 2]), A.ValidationError);
  assert.throws(() => M.statistics([1, NaN, 3], [1, 2, 3]), A.ValidationError);
  assert.throws(() => M.populationCorrelation(-1), A.ValidationError);
  console.log("Correlation: exact correlations, confounding limits, fixed populations, degenerate variance and validation passed.");
} finally { Math.random = originalRandom; }
