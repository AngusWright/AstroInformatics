"use strict";
const assert = require("node:assert/strict");
const A = require("../TeachingTools/statistics.js"), M = require("./model.js");
const close = (actual, expected, tolerance = 1e-11) => assert.ok(Math.abs(actual - expected) <= tolerance,
  `${actual} != ${expected} (tolerance ${tolerance})`);
const originalRandom = Math.random;
Math.random = () => { throw new Error("Use the seeded RNG, not Math.random"); };
try {
  // Retained model.R fit_colour_classifier: weight=c(4.5,3), threshold=2.25 for these fixed colours.
  const x = [[-2, -1], [-1, -2], [-1, 0], [0, -1], [1, 1], [2, 0], [2, 2], [3, 1]];
  const labels = [0, 0, 0, 0, 1, 1, 1, 1], fit = M.fitClassifier(x, labels);
  close(fit.weight[0], 4.5); close(fit.weight[1], 3); close(fit.threshold, 2.25);
  assert.deepEqual(M.predict(fit, [[0, 0], [1, 0], [0, 1], [0.5, 0]]), [0, 1, 1, 0],
    "Ties must follow the retained R rule and be assigned to class zero");
  assert.deepEqual(M.predict(fit, [[1.5, 0.6], [2.5, 0.6], [1.5, 1.6]], [1.5, 0.6]), [0, 1, 1]);
  const relabelled = M.fitClassifier(x, labels.map(label => 1 - label));
  close(relabelled.weight[0], -fit.weight[0]); close(relabelled.weight[1], -fit.weight[1]);
  close(relabelled.threshold, -fit.threshold);

  const base = M.prepare(42), matched = M.compute({ offset: 0, noise: 0 }, base);
  const shifted = M.compute({ offset: 1.5, noise: 0 }, base);
  assert.equal(matched.training, shifted.training); assert.equal(matched.model, shifted.model);
  assert.equal(matched.original, shifted.original); assert.equal(matched.originalPrediction, shifted.originalPrediction);
  close(matched.originalAccuracy, shifted.originalAccuracy);
  shifted.target.forEach((row, i) => { close(row[0] - matched.target[i][0], 1.5); close(row[1] - matched.target[i][1], 0.6); });
  assert.ok(matched.frozenAccuracy > 0.95 && shifted.frozenAccuracy < 0.8);
  assert.deepEqual(matched, M.compute({ offset: 0, noise: 0 }), "Seeded experiments must be repeatable");
  assert.notDeepEqual(base.training, M.prepare(43).training);
  const corrected = M.compute({ offset: 1.5, noise: 0, adaptation: "known" }, base);
  assert.deepEqual(corrected.adaptedPrediction, matched.frozenPrediction, "Known calibration must recover the matched-survey decisions");
  const noisy = M.compute({ offset: 1.5, noise: 1.2, adaptation: "known" }, base);
  assert.equal(noisy.model, base.model); assert.equal(noisy.originalPrediction, base.originalPrediction);
  assert.ok(noisy.adaptedAccuracy < corrected.adaptedAccuracy - 0.15, "A calibration correction must not silently remove noise");
  const aligned = M.compute({ offset: 1.5, noise: 0.3, adaptation: "unlabelled" }, base);
  const alignedOpposite = M.compute({ offset: -1.5, noise: 0.3, adaptation: "unlabelled" }, base);
  assert.deepEqual(aligned.adaptedPrediction, alignedOpposite.adaptedPrediction, "Mean alignment is equivariant to an additive offset");
  close(aligned.correction[0] - alignedOpposite.correction[0], 3);
  close(aligned.correction[1] - alignedOpposite.correction[1], 1.2);
  const hiddenTruth = { ...base, labels: base.labels.map(label => 1 - label) };
  const differentTruth = M.compute({ offset: 1.5, noise: 0.3, adaptation: "unlabelled" }, hiddenTruth);
  assert.deepEqual(differentTruth.correction, aligned.correction);
  assert.deepEqual(differentTruth.adaptedPrediction, aligned.adaptedPrediction, "Adaptation must not use target truth");
  close(differentTruth.adaptedAccuracy, 1 - aligned.adaptedAccuracy);
  assert.deepEqual(M.calibrationCorrection([[0, 0], [2, 4]], [[3, 2], [5, 6]], "unlabelled"), [3, 2]);
  const view = M.compute({ view: "original", offset: 1.5, noise: 0 }, base);
  assert.deepEqual(view.target, shifted.target); assert.equal(view.model, shifted.model);
  for (const seed of [0, 2147483647]) for (const offset of [-2, 2]) for (const noise of [0, 1.2]) {
    const result = M.compute({ seed, offset, noise, adaptation: "known" });
    assert.ok(result.target.flat().every(Number.isFinite));
    assert.ok([result.originalAccuracy, result.frozenAccuracy, result.adaptedAccuracy].every(value => value >= 0 && value <= 1));
  }
  for (const bad of [{ seed: -1 }, { seed: 0.5 }, { seed: 2147483648 }, { offset: 2.1 }, { offset: NaN },
    { noise: -0.1 }, { noise: Infinity }, { noise: 1.3 }, { view: "training" }, { adaptation: "target-labels" }]) {
    assert.throws(() => M.compute(bad), A.ValidationError);
  }
  assert.throws(() => M.fitClassifier([[0, 0], [1, 1], [2, 2], [3, 3]], [0, 0, 1, 1]), A.ValidationError);
  assert.throws(() => M.fitClassifier(x, labels.map(() => 0)), A.ValidationError);
  assert.throws(() => M.fitClassifier(x, [0, 1]), A.ValidationError);
  assert.throws(() => M.predict(fit, [[NaN, 0]]), A.ValidationError);
  assert.throws(() => M.calibrationCorrection(x, x, "labels"), A.ValidationError);
  console.log("SurveyShift: retained-R LDA, frozen source training, independent validation, honest corrections and validation passed.");
} finally { Math.random = originalRandom; }
