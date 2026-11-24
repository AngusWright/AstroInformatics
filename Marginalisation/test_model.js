"use strict";
const assert = require("node:assert/strict");
const A = require("../TeachingTools/statistics.js");
const M = require("./model.js");
const close = (actual, expected, tolerance = 1e-11) =>
  assert.ok(Math.abs(actual - expected) < tolerance, `${actual} != ${expected}`);

// Direct evaluations from retained model.R (R 17-digit output).
close(M.mixtureJoint(-1, 0, 0.65, 1), 0.91223017772492121);
close(M.mixtureMarginal(1.3, 0.65, 1), 0.25931248226093129);
close(M.mixtureInterval(0.9, 1.5, 0.65, 1), 0.15254294318212769);
close(M.mixtureInterval(0.9, 1.5, 0.7, 1.5), 0.11072054431171223);
close(M.mixtureInterval(-1.23, -0.77, 0, 1), 0.6826894921370859);
assert.equal(M.mixtureInterval(0, 0, 0.7, 1.5), 0);

// Separate midpoint integration checks the dy factor and true marginalisation.
const integrate = (fn, lo, hi, n) => {
  const dx = (hi - lo) / n;
  let total = 0;
  for (let i = 0; i < n; i++) total += fn(lo + (i + 0.5) * dx) * dx;
  return total;
};
for (const x of [-1, 0, 1.3, 3])
  close(integrate(y => M.mixtureJoint(x, y, 0.7, 1.5), -14, 16, 3000),
    M.mixtureMarginal(x, 0.7, 1.5), 1e-10);
close(integrate(x => M.mixtureMarginal(x, 0.7, 1.5), -14, 16, 6000), 1, 1e-10);
close(integrate(x => M.mixtureMarginal(x, 0.7, 1.5), -0.5, 0.8, 2000),
  M.mixtureInterval(-0.5, 0.8, 0.7, 1.5), 1e-8);
// Independent R integrate() over x with conditional-Normal y-window probabilities.
close(M.integrateGrid(0.7, 2).mass, 0.99759719104748223, 3e-6);
close(M.integrateGrid(0.7, 2, 300, 240).mass, 0.99759719104748223, 2e-6);
const p = { seed: 42, weight: 0.7, width: 1.5, centre: 1.2, show_profile: false };
const result = M.compute(p);
assert.ok(result.mode < 0, "The narrow mode is tallest despite the broad component holding most mass");
assert.ok(Math.abs(result.profile[60] - result.marginal[60]) > 0.01, "Profiling must not silently stand in for integration");
assert.deepEqual(result.joint, M.compute({ ...p, centre: -1, show_profile: true, seed: 43 }).joint);
assert.ok(M.compute({ ...p, width: 2, centre: 3, weight: 0.95 }).probability > 0);
for (const bad of [{ seed: 1.5 }, { weight: 1 }, { width: 0 }, { centre: Infinity }, { show_profile: "yes" }])
  assert.throws(() => M.compute({ ...p, ...bad }), A.ValidationError);
assert.throws(() => M.mixtureInterval(1, -1, 0.7, 1.5), A.ValidationError);
assert.throws(() => M.integrateGrid(0.5, 1, 0, 10), A.ValidationError);
console.log("Marginalisation: R densities, analytic probabilities, area-weighted integration and profiling distinction passed.");
