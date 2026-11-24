"use strict";
const assert = require("node:assert/strict");
const { test } = require("node:test");
const A = require("../TeachingTools/statistics.js");
const M = require("./model.js");
const base = { seed: 42, position: 8, contaminants: 1, n: 40, trim: 0.1 };
const close = (actual, expected, tolerance = 1e-12) =>
  assert.ok(Math.abs(actual - expected) <= tolerance, `${actual} != ${expected} (tolerance ${tolerance})`);

test("mean, even/odd median, sample SD and scaled MAD match independent R answers", () => {
  for (const [values, expected] of [
    [[-2, -1, 0, 1, 100], { mean: 19.6, median: 0, sd: 44.958870092563494, nmad: 1.482602218505602 }],
    [[1, 2, 3, 100], { mean: 26.5, median: 2.5, sd: 49.00680224893955, nmad: 1.482602218505602 }]
  ]) {
    const summaries = M.summaries(values);
    for (const key of Object.keys(expected)) close(summaries[key], expected[key], 2e-12);
  }
  const tied = M.summaries([3, 3, 3, 3]);
  assert.equal(tied.mean, 3); assert.equal(tied.median, 3);
  assert.equal(tied.sd, 0); assert.equal(tied.nmad, 0);
});

test("trimmed means remove floor(n*t) values from each tail, as R does", () => {
  const values = [-10, -2, 0, 1, 3, 8, 20], original = [...values];
  close(M.trimmedMean(values, 0), 20 / 7);
  close(M.trimmedMean(values, 0.1), 20 / 7);
  close(M.trimmedMean(values, 0.15), 2);
  close(M.trimmedMean(values, 0.3), 4 / 3);
  close(M.trimmedMean(values, 0.49), 1);
  assert.deepEqual(values, original, "Computing a trimmed mean must not reorder observations");
});

test("an isolated extreme changes mean and SD but not median, MAD or a suitable trim", () => {
  const moderate = M.summaries([-2, -1, 0, 1, 100], 0.2);
  const extreme = M.summaries([-2, -1, 0, 1, 1000], 0.2);
  close(extreme.mean - moderate.mean, 180);
  assert.ok(extreme.sd > 9 * moderate.sd);
  assert.equal(moderate.median, extreme.median);
  assert.equal(moderate.nmad, extreme.nmad);
  assert.equal(moderate.trimmed, 0); assert.equal(extreme.trimmed, 0);
});

test("contamination replaces values, preserves offsets, and has analytic mean response c/n", () => {
  const clean = [-2, -1, 0, 1, 2], offsets = [0.1, -0.2];
  assert.deepEqual(M.contaminate(clean, offsets, 2, 10), [10.1, 9.8, 0, 1, 2]);
  assert.deepEqual(M.contaminate(clean, offsets, 0, 10), clean);
  assert.deepEqual(clean, [-2, -1, 0, 1, 2]);
  const data = M.prepare(42);
  for (const contaminants of [0, 1, 10]) {
    const result = M.compute({ ...base, contaminants }, data);
    close(result.influence.at(-1).mean - result.influence[0].mean, 24 * contaminants / base.n);
    if (contaminants === 0) {
      assert.deepEqual(result.observed, result.clean);
      assert.ok(result.influence.every(row => row.mean === result.summaries.mean));
    }
  }
  const half = M.compute({ ...base, n: 20, contaminants: 10, position: 12 }, data);
  assert.ok(half.summaries.median > 4, "A median is not immune to half-sample contamination");
});

test("moving outliers, trimming and sample-size changes do not redraw clean observations", () => {
  const data = M.prepare(42), first = M.compute(base, data);
  const moved = M.compute({ ...base, position: -12, trim: 0.4 }, data);
  const larger = M.compute({ ...base, n: 120, contaminants: 10 }, data);
  assert.strictEqual(first.data, moved.data);
  assert.deepEqual(first.clean, moved.clean);
  assert.deepEqual(first.jitter, moved.jitter);
  assert.deepEqual(first.observed.slice(1), moved.observed.slice(1));
  close(first.observed[0] - moved.observed[0], 20);
  assert.deepEqual(first.clean, larger.clean.slice(0, 40));
  assert.deepEqual(first.jitter, larger.jitter.slice(0, 40));
  assert.deepEqual(M.prepare(42), data);
  assert.notDeepEqual(M.prepare(43).clean, data.clean);
});

test("invalid data, trims and controls produce validation errors", () => {
  for (const patch of [{ seed: -1 }, { seed: 0.5 }, { n: 19 }, { n: 121 }, { n: 20.5 },
    { contaminants: -1 }, { contaminants: 11 }, { contaminants: 0.5 }, { position: -12.1 },
    { position: Infinity }, { trim: -0.1 }, { trim: 0.41 }]) {
    assert.throws(() => M.compute({ ...base, ...patch }), A.ValidationError);
  }
  for (const values of [[], [1], [1, NaN]]) assert.throws(() => M.summaries(values), A.ValidationError);
  for (const trim of [-0.1, 0.5, NaN]) assert.throws(() => M.trimmedMean([1, 2], trim), A.ValidationError);
  assert.throws(() => M.contaminate([1], [], 1, 0), A.ValidationError);
  assert.throws(() => M.compute(base, M.prepare(0)), A.ValidationError);
});
