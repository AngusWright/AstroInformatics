"use strict";
const test = require("node:test");
const assert = require("node:assert/strict");
const Astro = require("../TeachingTools/statistics.js");
const M = require("./model.js");
const config = overrides => ({
  family: "gaussian", noise: "log-linear", seed: 666, samples: 81, repeats: 12,
  xmin: -4, xmax: 4, noiseMin: -2, noiseMax: -1, ...overrides
});
const settings = overrides => ({
  widthMin: .9, widthMax: 1.1, models: 3, weighting: "legacy", residual: false, ...overrides
});
const close = (actual, expected, tolerance = 1e-9) =>
  assert.ok(Math.abs(actual - expected) <= tolerance, `${actual} != ${expected} (tolerance ${tolerance})`);

test("seeded observations are repeatable, use the actual noise SD, and do not consume global randomness", () => {
  const original = Math.random;
  Math.random = () => { throw Error("Global RNG used"); };
  try {
    const a = M.createData(config()), b = M.createData(config());
    assert.deepEqual(a, b);
    assert.notDeepEqual(a.observations, M.createData(config({ seed: 667 })).observations);
    close(a.sigma[0], .01); close(a.sigma.at(-1), .1);
    close(a.sigma[40], Math.sqrt(.001));
    const snapshot = JSON.stringify(a);
    M.compare(a, settings());
    M.compare(a, settings({ weighting: "known", widthMin: .8 }));
    assert.equal(JSON.stringify(a), snapshot, "Inference mutated or resimulated observations");
  } finally { Math.random = original; }
});

test("all noise-generation branches preserve the R endpoint definitions", () => {
  const a = M.createData(config({ noise: "constant" }));
  a.sigma.forEach(s => close(s, .055));
  const b = M.createData(config({ noise: "linear" }));
  close(b.sigma[40], .055);
  const c = M.createData(config({ noiseMin: -1, noiseMax: -1 }));
  c.sigma.forEach(s => close(s, .1));
  const d = M.createData(config({ family: "linear", samples: 3, xmin: -1, xmax: 1 }));
  assert.deepEqual(d.truth, [-1, 0, 1]);
});

test("linear intercept known answers and exact heteroskedastic propagated SD", () => {
  const x = [-1, 0, 1], y = [0, 2, 4], sigma = [1, 2, 4];
  const w = M.weights(sigma, "known");
  const fitted = M.fit(x, y, sigma, 2, "linear", w);
  close(fitted.shift, 2);
  close(fitted.se, 1 / Math.sqrt(1 + .25 + .0625));
  const equal = M.fit(x, y, sigma, 2, "linear", M.weights(sigma, "equal"));
  close(equal.se, Math.sqrt(21) / 3);
  const legacy = M.fit(x, y, sigma, 2, "linear", M.weights(sigma, "legacy"));
  close(legacy.se, Math.sqrt(3) / 1.75);
});

test("wrong slope creates a weight-dependent noiseless bias, not a rescaling of uncertainty", () => {
  const x = [-1, 0, 1], sigma = [.1, .2, .4];
  const known = M.fit(x, x, sigma, .8, "linear", M.weights(sigma, "known"));
  const legacy = M.fit(x, x, sigma, .8, "linear", M.weights(sigma, "legacy"));
  const equal = M.fit(x, x, sigma, .8, "linear", M.weights(sigma, "equal"));
  close(known.shift, .2 * (-100 + 6.25) / 131.25);
  close(legacy.shift, .2 * (-10 + 2.5) / 17.5);
  close(equal.shift, 0);
  assert.ok(known.se < legacy.se && legacy.se < equal.se);
});

test("Gaussian amplitude profiling recovers a known shifted, scaled density", () => {
  const x = M.grid(-4, 4, 121), sigma = x.map(() => .03);
  const y = x.map(t => 1.7 * Astro.normalPdf(t, .35, .8));
  const w = M.weights(sigma, "known"), fitted = M.fit(x, y, sigma, .8, "gaussian", w);
  close(fitted.shift, .35, 2e-7);
  close(fitted.amplitude, 1.7, 2e-7);
  assert.equal(fitted.boundary, false);
  let a = 0, b = 0, c = 0;
  for (let i = 0; i < x.length; i++) {
    const g = Astro.normalPdf(x[i], .35, .8), d = 1.7 * g * (x[i] - .35) / (.8 ** 2);
    a += w[i] * d * d; b += w[i] * d * g; c += w[i] * g * g;
  }
  close(fitted.se, Math.sqrt(c / (a * c - b * b)), 1e-8);
  close(fitted.amplitudeSE, Math.sqrt(a / (a * c - b * b)), 1e-8);
});

test("misspecified Gaussian uncertainty includes residual curvature: finite-difference sensitivity", () => {
  const d = M.createData(config({ samples: 25, repeats: 1 }));
  const w = M.weights(d.sigma, "known"), width = .75;
  const f = M.fit(d.x, d.truth, d.sigma, width, "gaussian", w);
  assert.ok(Math.abs(f.shift) > .01, "Asymmetric weighting should move the wrong-width Gaussian");
  let variance = 0;
  const eps = .0003;
  for (let i = 0; i < d.x.length; i++) {
    const yp = d.truth.slice(), ym = d.truth.slice(); yp[i] += eps; ym[i] -= eps;
    const plus = M.fit(d.x, yp, d.sigma, width, "gaussian", w).shift;
    const minus = M.fit(d.x, ym, d.sigma, width, "gaussian", w).shift;
    variance += (((plus - minus) / (2 * eps)) * d.sigma[i]) ** 2;
  }
  close(f.se, Math.sqrt(variance), .00005);
});

test("Gaussian support boundary and singular fits never receive symmetric confidence intervals", () => {
  const x = M.grid(-1, 1, 51), sigma = x.map(() => .01), w = M.weights(sigma, "known");
  const f = M.fit(x, x.map(t => Astro.normalPdf(t, 2, .7)), sigma, .7, "gaussian", w);
  assert.equal(f.boundary, true);
  assert.equal(f.se, null);
  const flat = M.fit(x, x.map(() => 0), sigma, 1, "gaussian", w);
  assert.equal(flat.se, null);
});

test("linear normal intervals cover the reference, not necessarily the misspecified truth", () => {
  const d = M.createData(config({ family: "linear", samples: 17, repeats: 1600,
    xmin: 0, xmax: 4, noiseMin: -1.5, noiseMax: -1 }));
  const groups = M.compare(d, settings({ widthMin: .5, widthMax: 1, models: 2, weighting: "known" }));
  for (const g of groups) {
    assert.equal(g.summary.intervalCount, 1600);
    assert.ok(g.summary.coverageReference > .925 && g.summary.coverageReference < .975);
    assert.ok(Math.abs(g.summary.sd / g.summary.referenceSE - 1) < .07);
  }
  close(groups[0].summary.coverageTruth, 0);
  assert.ok(groups[1].summary.coverageTruth > .925);
});

test("residual heuristic fixes the linear slope call, reuses data, and makes no analytic coverage claim", () => {
  const d = M.createData(config({ family: "linear", noise: "constant", samples: 61, repeats: 8 }));
  const before = JSON.stringify(d);
  const groups = M.compare(d, settings({ residual: true }));
  assert.equal(JSON.stringify(d), before);
  assert.ok(groups.every(g => g.fits.every(f => f.se === null && Number.isFinite(f.shift))));
  assert.ok(groups.every(g => g.summary.coverageTruth === null && g.summary.intervalCount === 0));
  assert.ok(groups.every(g => g.reference.se > 0));
  const x = [0, 1, 2, 3], y = [0, 0, 0, 10], base = x.map(() => 1);
  assert.throws(() => M.residualWeights(x, y, { shift: 0, amplitude: null }, 1, "linear", base), Astro.ValidationError);
});

test("Gaussian residual correction profiles amplitude and reports empirical-only noisy-fit uncertainty", () => {
  const d = M.createData(config({ noise: "constant", repeats: 8 }));
  const groups = M.compare(d, settings({ residual: true }));
  assert.ok(groups.every(g => g.fits.every(f => Number.isFinite(f.amplitude) && f.amplitudeSE === null && f.se === null)));
  assert.ok(groups.every(g => g.summary.sd > 0 && g.summary.coverageReference === null));
});

test("paired trend uncertainty uses independent realisations, not model count", () => {
  const groups = [1, 2, 3].map(width => ({ width, fits: [0, 1, 2, 3].map(z => ({ shift: 2 * width + z })) }));
  const trend = M.biasTrend(groups);
  close(trend[0].mean, 3.5);
  close(trend.at(-1).mean, 7.5);
  close(trend[0].sd, Math.sqrt(5 / 3));
  close(trend[0].meanSE, Math.sqrt(5 / 3) / 2);
  assert.throws(() => M.biasTrend(groups.map(g => ({ ...g, width: 1 }))), Astro.ValidationError);
});

test("one realisation and repeated model widths have honest unavailable empirical summaries", () => {
  const d = M.createData(config({ family: "linear", repeats: 1 }));
  const groups = M.compare(d, settings({ widthMin: 1, widthMax: 1 }));
  assert.equal(groups[0].summary.sd, null);
  assert.equal(groups[0].summary.mcse, null);
  assert.deepEqual(groups[0].fits, groups[2].fits);
  const h = M.histogram(groups[0].fits.map(f => f.shift));
  close(h.points.reduce((s, p) => s + p[1] * h.step, 0), 1);
});

test("unsupported, reversed, nonfinite, underidentified and oversized inputs fail explicitly", () => {
  for (const change of [
    { family: "gamma" }, { family: "poissonian" }, { noise: "oops" }, { seed: -1 },
    { seed: 1.2 }, { samples: 1 }, { samples: 2 }, { repeats: 0 }, { xmin: 4, xmax: -4 },
    { noiseMin: -1, noiseMax: -2 }, { noiseMax: NaN }, { xmin: -11 }, { samples: 10000, repeats: 10000 }
  ]) assert.throws(() => M.createData(config(change)), Astro.ValidationError);
  const d = M.createData(config());
  for (const change of [{ models: 1 }, { widthMin: 2 }, { widthMin: 1.2, widthMax: 1 }, { weighting: "wrong" }])
    assert.throws(() => M.compare(d, settings(change)), Astro.ValidationError);
  const large = M.createData(config({ samples: 1000, repeats: 400 }));
  assert.throws(() => M.compare(large, settings()), /budget/);
  assert.throws(() => M.fit([0, 1], [1], [1, 1], 1, "linear", [1, 1]), Astro.ValidationError);
  assert.throws(() => M.fit([0, 1], [1, 2], [1, 0], 1, "linear", [1, 1]), Astro.ValidationError);
});
