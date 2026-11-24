"use strict";
const assert = require("node:assert/strict");
const { test } = require("node:test");
const M = require("./model.js");
const { ValidationError } = require("../TeachingTools/statistics.js");

test("magnitude selection matches the R relation, including the detection edge", () => {
  assert.equal(M.observed(0, 100, 17), true);
  assert.equal(M.observed(0, 100, 16.999), false);
  assert.equal(M.threshold(100, 20), -1.2);
  // R: selection_observed(c(-2,-1.2,0), rep(100,3), 20).
  assert.deepEqual([-2, -1.2, 0].map(l => M.observed(l, 100, 20)), [false, true, true]);
  assert.ok(Math.abs(M.threshold(100, 20) - M.threshold(10, 20) - 2) < 1e-14);
});

test("selection and reveal use one reproducible parent population", () => {
  const shallow = M.compute({ limit: 17 }), deep = M.compute({ limit: 24 });
  assert.deepEqual(shallow, M.compute({ limit: 17 }));
  assert.deepEqual(shallow.distance, deep.distance);
  assert.deepEqual(shallow.luminosity, deep.luminosity);
  assert.deepEqual(deep, M.compute({ limit: 24, overlay: "hide" }));
  assert.ok(shallow.distance.every(d => d >= 30 && d <= 800));
  assert.ok(shallow.final.count < deep.final.count);
  shallow.final.detected.forEach((yes, i) => { if (yes) assert.equal(deep.final.detected[i], true); });
  for (const phase of [0, 0.5, 1]) {
    const state = M.atPhase(shallow, phase);
    assert.equal(state.count + state.missing, 2400);
    assert.equal(state.binCounts.reduce((a, b) => a + b, 0), state.count);
    assert.equal(state.limit, 26 + phase * (17 - 26));
  }
  assert.deepEqual(M.atPhase(shallow), shallow.final);
  assert.notDeepEqual(shallow.distance, M.compute({ seed: 43 }).distance);
});

test("no detections produce no invented means; empty bins remain empty", () => {
  const empty = M.atPhase({ distance: [40, 120], luminosity: [-10, -10], targetLimit: 17 });
  assert.equal(empty.count, 0);
  assert.equal(empty.mean, null);
  assert.ok(empty.binMeans.every(value => value === null));
  const one = M.atPhase({ distance: [80, 160], luminosity: [0, -10], targetLimit: 17 });
  assert.equal(one.binCounts[0], 1);
  assert.equal(one.binMeans[0], 0);
  assert.equal(one.binMeans[1], null);
});

test("pure API rejects invalid domains with expected validation errors", () => {
  for (const p of [{ limit: 16 }, { limit: NaN }, { overlay: "unknown" }, { seed: -1 }, { seed: 1.2 }]) {
    assert.throws(() => M.compute(p), ValidationError);
  }
  assert.throws(() => M.observed(0, 0, 20), ValidationError);
  assert.throws(() => M.atPhase(M.compute(), 2), ValidationError);
  assert.doesNotThrow(() => M.compute({ seed: 2147483647 }));
});
