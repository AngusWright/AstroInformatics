"use strict";
const { test } = require("node:test");
const assert = require("node:assert/strict");
const A = require("./statistics.js");
function near(actual, expected, tolerance) {
  assert.ok(Math.abs(actual - expected) <= tolerance, `${actual} != ${expected} (+/- ${tolerance})`);
}

test("Beta and Gamma quantiles match independent R references, including boundary concentration", () => {
  const probabilities = [0.001, 0.025, 0.5, 0.975, 0.999];
  // qbeta/qgamma from the retained R implementation, not inverse-CDF self checks.
  const cases = [
    [A.beta(0.5, 50), [1.5786706587275014e-8, 9.8698680204510812e-6, 0.0045617224731038649, 0.049236196384272457, 0.10310226341871376]],
    [A.beta(9, 8), [0.1899303961785927, 0.29877689908317751, 0.53059173132260329, 0.75348988850942467, 0.85126011791932188]],
    [A.beta(10000.5, 0.5), [0.99945878173783509, 0.99974884351561666, 0.99997725400620252, 0.99999995089777294, 0.99999999992146216]],
    [A.gamma(0.5, 0.1), [7.8539857463124513e-6, 0.0049103455858762806, 2.274682115597864, 25.119430936574421, 54.137830853313652]],
    [A.gamma(43, 12), [2.1285250918822047, 2.5932761009972647, 3.555594162399903, 4.7309832315408871, 5.5115551056227456]],
    [A.gamma(100050, 100020), [0.99055576501591303, 0.99411115511112624, 1.0002966073471715, 1.0065076641702089, 1.0101011004356382]]
  ];
  for (const [distribution, expected] of cases) probabilities.forEach((p, i) => {
    const q = distribution.quantile(p);
    near(q, expected[i], Math.abs(expected[i]) * 1e-8 + 1e-20);
    near(distribution.cdf(q), p, 2e-8);
  });
});

test("Normal tails and quantiles retain precision", () => {
  const xs = [-8, -2, 0, 1.96, 8];
  const cdfs = [6.2209605742717849e-16, 0.022750131948179212, 0.5, 0.97500210485177963, 0.99999999999999933];
  xs.forEach((x, i) => near(A.normalCdf(x), cdfs[i], Math.abs(cdfs[i]) * 1e-12 + 1e-30));
  near(A.normalQuantile(0.025), -1.959963984540054, 1e-12);
  near(A.normalQuantile(0.975, 4, 2), 7.919927969080108, 1e-12);
  near(A.normalCdf(-Infinity), 0, 0); near(A.normalCdf(Infinity), 1, 0);
});

test("distributions respect support and proper parameterisation", () => {
  assert.equal(A.gamma(1, 2).density(0), 2);
  assert.equal(A.gamma(0.5, 2).density(0), Infinity);
  assert.equal(A.beta(1, 3).density(0), 3);
  assert.equal(A.beta(3, 1).density(1), 3);
  assert.equal(A.beta(0.5, 0.5).density(0), Infinity);
  assert.equal(A.gamma(3, 2).density(-1), 0);
  near(A.gamma(3, 2).mean, 1.5, 0);
  near(A.gamma(3, 2).sd ** 2, 0.75, 1e-15);
  for (const action of [() => A.beta(0, 1), () => A.gamma(1, -2), () => A.normal(0, 0),
    () => A.beta(1, 1).quantile(-0.1), () => A.gamma(1, 1).quantile(1.1)]) {
    assert.throws(action, A.ValidationError);
  }
});

test("sample quantiles interpolate like R type 7 without reordering observations", () => {
  const observations = [10, 1, 5, 2];
  for (const [probability, expected] of [[0, 1], [0.1, 1.3], [0.25, 1.75], [0.5, 3.5], [0.9, 8.5], [1, 10]]) {
    near(A.quantile(observations, probability), expected, 1e-14);
  }
  assert.deepEqual(observations, [10, 1, 5, 2]);
  assert.throws(() => A.quantile([], 0.5), A.ValidationError);
});

test("histograms close the last bin and retain omitted tail probability", () => {
  const histogram = A.histogram([-2, 0, 1, 1.5, 2, 3, 4], [0, 1, 2, 4]);
  assert.deepEqual(histogram.counts, [1, 2, 3]);
  assert.equal(histogram.excluded, 1);
  near(histogram.density.reduce((sum, height, i) => sum + height * histogram.widths[i], 0), 6 / 7, 1e-15);
  assert.deepEqual(A.histogram([], [0, 1, 2]).density, [0, 0]);
  assert.throws(() => A.histogram([1], [0, 0, 2]), A.ValidationError);
});
