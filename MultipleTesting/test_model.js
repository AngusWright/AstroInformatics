"use strict";
const assert = require("node:assert/strict");
const A = require("../TeachingTools/statistics.js");
const M = require("./model.js");
const close = (actual, expected, tolerance = 1e-12) =>
  assert.ok(Math.abs(actual - expected) <= tolerance * Math.max(Math.abs(expected), 1e-300), `${actual} != ${expected}`);

// Retained model.R / pnorm references, independent of the browser generator.
const reference = M.searchFalseProbability(3, 1600);
close(reference.single, 0.0013498980316300946);
close(reference.global, 0.88482423701342394);
close(reference.expected, 2.1598368506081513);
close(M.searchFalseProbability(5, 3600, true).global, 0.002061763548728873);
close(M.searchFalseProbability(3, 1).global, reference.single);
assert.equal(M.searchFalseProbability(3, 0).global, 0);
assert.equal(M.searchFalseProbability(0, 0, true).global, 0);
assert.equal(M.searchFalseProbability(0, 1, true).global, 1);
// A subtraction 1 - Phi(10) would round to zero.
close(M.searchFalseProbability(10, 1000000000).global, 7.619853024160495e-15, 2e-12);

const handSurvey = { threshold: 3, twoSided: false, cells: [
  { truth: false, hit: true }, { truth: true, hit: true },
  { truth: false, hit: false }, { truth: true, hit: false }
] };
const hand = M.reveal(handSurvey);
assert.deepEqual([hand.trueDiscoveries, hand.falseDiscoveries, hand.nullTrials], [1, 1, 2]);
assert.equal(hand.falseDiscoveryFraction, 0.5);
assert.equal(M.reveal(handSurvey, 0).falseDiscoveryFraction, null);
assert.equal(M.reveal(handSurvey, 50).tested, 2);

const p = { seed: 42, side: 40, threshold: 3, tails: "one", signal_fraction: 0, amplitude: 4 };
const originalRandom = Math.random;
Math.random = () => { throw new Error("Unseeded randomness is forbidden"); };
try {
  const first = M.compute(p), lowered = M.compute({ ...p, threshold: 1.5, tails: "two" });
  assert.deepEqual(first.cells.map(c => c.z), lowered.cells.map(c => c.z));
  assert.deepEqual(first.cells, M.compute(p).cells);
  assert.deepEqual(first.cells, M.compute({ ...p, side: 60 }).cells.slice(0, first.cells.length));
  assert.notDeepEqual(first.cells, M.compute({ ...p, seed: 43 }).cells);
  assert.equal(M.reveal(first).trueDiscoveries, 0);
  assert.ok(M.reveal(lowered).falseDiscoveries > 0);
  const injected = M.reveal(M.compute({ ...p, signal_fraction: 0.1, amplitude: 8 }));
  assert.ok(injected.trueDiscoveries > 0);
  assert.equal(injected.nullTrials + injected.trueDiscoveries, 1600);
  M.compute({ ...p, seed: 2147483647, side: 5, threshold: 5 });
} finally { Math.random = originalRandom; }

// Independent repeated-survey calibration, also used in test_numerics.R.
const rng = A.createRng(701), surveys = 20000;
let falseAlarms = 0;
for (let run = 0; run < surveys; run++) {
  let hit = false;
  for (let i = 0; i < 20; i++) if (rng.normal() > 2.5) hit = true;
  if (hit) falseAlarms++;
}
assert.ok(Math.abs(falseAlarms / surveys - M.searchFalseProbability(2.5, 20).global) < 0.01);
for (const bad of [{ seed: -1 }, { side: 4 }, { threshold: NaN }, { tails: "invalid" }, { signal_fraction: -0.1 }])
  assert.throws(() => M.compute({ ...p, ...bad }), A.ValidationError);
assert.throws(() => M.searchFalseProbability(3, 1.5), A.ValidationError);
assert.throws(() => M.reveal(handSurvey, 101), A.ValidationError);
console.log("MultipleTesting: R references, stable tails, FWER calibration, discovery accounting and fixed-sky controls passed.");
