(function (root, factory) {
  "use strict";
  const api = factory(typeof module === "object" && module.exports ? require("../TeachingTools/statistics.js") : root.Astro);
  if (typeof module === "object" && module.exports) module.exports = api;
  root.SurveyShiftModel = api;
})(globalThis, function (A) {
  "use strict";
  const defaults = { seed: 42, offset: 0.8, noise: 0.3, view: "target", adaptation: "none" }, SIZE = 400;
  function validate(p) {
    A.check(Number.isInteger(p.seed) && p.seed >= 0 && p.seed <= 2147483647, "Enter an integer seed between 0 and 2147483647.");
    A.check(Number.isFinite(p.offset) && p.offset >= -2 && p.offset <= 2, "Calibration offset must be between -2 and 2.");
    A.check(Number.isFinite(p.noise) && p.noise >= 0 && p.noise <= 1.2, "Additional noise must be between 0 and 1.2.");
    A.check(["target", "original"].includes(p.view), "Choose the original or new validation survey.");
    A.check(["none", "known", "unlabelled"].includes(p.adaptation), "Choose a valid calibration correction.");
  }
  function mean2(x) {
    A.check(x.length > 0 && x.every(row => row.length === 2 && row.every(Number.isFinite)), "Provide a nonempty two-colour sample.");
    return [0, 1].map(j => x.reduce((sum, row) => sum + row[j] / x.length, 0));
  }
  function fitClassifier(x, labels) {
    mean2(x);
    A.check(labels.length === x.length && labels.every(label => label === 0 || label === 1), "Training labels must be zero or one and match the sample.");
    const groups = [0, 1].map(label => x.filter((_, i) => labels[i] === label));
    A.check(groups.every(group => group.length >= 2), "Each training class needs at least two observations.");
    const means = groups.map(mean2);
    const covariances = groups.map((group, i) => {
      let xx = 0, yy = 0, xy = 0;
      group.forEach(row => {
        const dx = row[0] - means[i][0], dy = row[1] - means[i][1];
        xx += dx * dx; yy += dy * dy; xy += dx * dy;
      });
      return [xx / (group.length - 1), xy / (group.length - 1), yy / (group.length - 1)];
    });
    const [a, b, d] = covariances[0].map((value, j) => (value + covariances[1][j]) / 2);
    const determinant = a * d - b * b;
    A.check(Number.isFinite(determinant) && a > 0 && d > 0 && determinant > Number.EPSILON * a * d,
      "The pooled colour covariance is singular; a linear discriminant is not identifiable.");
    const dx = means[1][0] - means[0][0], dy = means[1][1] - means[0][1];
    const weight = [(d * dx - b * dy) / determinant, (a * dy - b * dx) / determinant];
    const threshold = weight.reduce((sum, value, j) => sum + value * (means[0][j] + means[1][j]) / 2, 0);
    A.check([...weight, threshold].every(Number.isFinite), "The classifier exceeds numerical precision.");
    return { weight, threshold, means, covariance: [[a, b], [b, d]] };
  }
  function predict(model, x, correction = [0, 0]) {
    A.check(model.weight.length === 2 && [...model.weight, model.threshold, ...correction].every(Number.isFinite) && correction.length === 2,
      "The classifier and correction must be finite.");
    A.check(x.every(row => row.length === 2 && row.every(Number.isFinite)), "Validation colours must be finite pairs.");
    return x.map(row => Number((row[0] - correction[0]) * model.weight[0] +
      (row[1] - correction[1]) * model.weight[1] > model.threshold));
  }
  function calibrationCorrection(training, target, mode, knownOffset) {
    if (mode === "none") return [0, 0];
    if (mode === "known") {
      A.check(Number.isFinite(knownOffset), "Known calibration offset must be finite.");
      return [knownOffset, 0.4 * knownOffset];
    }
    A.check(mode === "unlabelled", "Choose a valid calibration correction.");
    const sourceMean = mean2(training), targetMean = mean2(target);
    return targetMean.map((value, j) => value - sourceMean[j]);
  }
  function accuracy(prediction, labels) {
    A.check(prediction.length === labels.length && labels.length > 0, "Validation predictions and truth must match.");
    return prediction.filter((value, i) => value === labels[i]).length / labels.length;
  }
  function prepare(seed) {
    const rng = A.createRng(seed), labels = Array.from({ length: SIZE }, (_, i) => Number(i >= SIZE / 2));
    function sample() {
      const x = labels.map(label => (2 * label - 1) * 0.8 + 0.35 * rng.normal());
      return labels.map((label, i) => [x[i], (2 * label - 1) * 0.4 + 0.45 * rng.normal()]);
    }
    const training = sample(), original = sample(), targetIntrinsic = sample();
    const noiseX = labels.map(() => rng.normal()), noise = labels.map((_, i) => [noiseX[i], rng.normal()]);
    const model = fitClassifier(training, labels), originalPrediction = predict(model, original);
    return { seed, labels, training, original, targetIntrinsic, noise, model, originalPrediction,
      originalAccuracy: accuracy(originalPrediction, labels) };
  }
  function compute(parameters, prepared) {
    const p = { ...defaults, ...parameters };
    validate(p);
    const base = prepared || prepare(p.seed);
    A.check(base.seed === p.seed, "The cached training sample must match the seed.");
    const target = base.targetIntrinsic.map((row, i) => row.map((value, j) =>
      value + (j === 0 ? p.offset : 0.4 * p.offset) + p.noise * base.noise[i][j]));
    const correction = calibrationCorrection(base.training, target, p.adaptation, p.offset);
    const frozenPrediction = predict(base.model, target), adaptedPrediction = predict(base.model, target, correction);
    return { base, ...base, target, correction, frozenPrediction, adaptedPrediction,
      frozenAccuracy: accuracy(frozenPrediction, base.labels), adaptedAccuracy: accuracy(adaptedPrediction, base.labels) };
  }
  return { defaults, SIZE, validate, fitClassifier, predict, calibrationCorrection, accuracy, prepare, compute };
});
