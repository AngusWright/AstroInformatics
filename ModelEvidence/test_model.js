"use strict";
const assert = require("node:assert/strict");
const { test } = require("node:test");
const M = require("./model.js");
const { ValidationError } = require("../TeachingTools/statistics.js");
const close = (actual, expected, tolerance = 2e-12) => assert.ok(Math.abs(actual - expected) < tolerance,
  `${actual} != ${expected}`);

test("normalised Bayes factor agrees with R and an exact null-estimate limit", () => {
  const r = M.lineEvidence(0.7, 0.4, 2);
  close(r.logBf, -0.15669249977997257);
  close(r.probability, 0.46090682872892397);
  close(r.logMaximumRatio, 1.53125);
  const zero = M.lineEvidence(0, 1, Math.sqrt(3));
  close(Math.exp(zero.logBf), 0.5);
  close(zero.probability, 1 / 3);
  close(zero.logMaximumRatio, 0);
  close(M.lineEvidence(-0.7, 0.4, 2).logBf, r.logBf);
  close(M.lineEvidence(0.7, 0.4, 1e-7).logBf, 0, 1e-12);
});

test("evidence integrates likelihood times a unit-area proper prior, not best fits", () => {
  const estimate = 0.7, error = 0.4, prior = 2, e = M.lineEvidence(estimate, error, prior);
  const normal = (x, mean, sd) => Math.exp(-0.5 * ((x - mean) / sd) ** 2) / (sd * Math.sqrt(2 * Math.PI));
  const n = 16000, lo = -16, hi = 16, h = (hi - lo) / n;
  const fn = a => normal(estimate, a, error) * normal(a, 0, prior);
  let total = fn(lo) + fn(hi);
  for (let i = 1; i < n; i++) total += (i % 2 ? 4 : 2) * fn(lo + i * h);
  close(total * h / 3, Math.exp(e.logAlternative));
  close(Math.log(total * h / 3 / normal(estimate, 0, error)), e.logBf);
  assert.ok(e.logBf < e.logMaximumRatio);
});

test("full-spectrum evidence has the correct Gaussian normalisation", () => {
  const r = M.compute(), p = M.defaults;
  const directNull = r.flux.reduce((sum, flux) => sum - Math.log(p.noise * Math.sqrt(2 * Math.PI)) -
    (flux - 1) ** 2 / (2 * p.noise ** 2), 0);
  close(r.logEvidenceNull, directNull, 1e-11);
  close(r.logEvidenceLine - r.logEvidenceNull, r.evidence.logBf, 1e-11);
  // Independent multivariate-normal calculation via the determinant lemma:
  // covariance = sigma² I + w² g gᵀ.
  const gg = r.template.reduce((sum, g) => sum + g * g, 0);
  const gy = r.template.reduce((sum, g, i) => sum + g * (r.flux[i] - 1), 0);
  const yy = r.flux.reduce((sum, flux) => sum + (flux - 1) ** 2, 0);
  const logDet = r.flux.length * Math.log(p.noise ** 2) + Math.log1p(p.prior ** 2 * gg / p.noise ** 2);
  const quadratic = yy / p.noise ** 2 - p.prior ** 2 * gy ** 2 /
    (p.noise ** 2 * (p.noise ** 2 + p.prior ** 2 * gg));
  close(r.logEvidenceLine, -0.5 * (r.flux.length * Math.log(2 * Math.PI) + logDet + quadratic), 1e-11);
});

test("changing only the prior keeps observations and best fit fixed", () => {
  const narrow = M.compute({ prior: 0.1 }), broad = M.compute({ prior: 10 });
  assert.deepEqual(narrow, M.compute({ prior: 0.1 }));
  assert.deepEqual(narrow.flux, broad.flux);
  assert.equal(narrow.estimate, broad.estimate);
  assert.equal(narrow.evidence.logMaximumRatio, broad.evidence.logMaximumRatio);
  assert.notEqual(narrow.evidence.logBf, broad.evidence.logBf);
  assert.notDeepEqual(narrow.flux, M.compute({ seed: 43 }).flux);
  for (const amplitude of [-2, 3]) for (const noise of [0.3, 2]) for (const prior of [0.1, 10]) {
    const r = M.compute({ amplitude, noise, prior });
    assert.ok(Number.isFinite(r.evidence.logBf));
    assert.ok(r.evidence.probability >= 0 && r.evidence.probability <= 1);
    assert.ok(r.evidence.logBf <= r.evidence.logMaximumRatio);
  }
});

test("pure calls reject improper priors and invalid data settings", () => {
  for (const p of [{ amplitude: 4 }, { noise: 0 }, { prior: 0 }, { prior: Infinity }, { seed: -1 }]) {
    assert.throws(() => M.compute(p), ValidationError);
  }
  assert.throws(() => M.lineEvidence(1, 0, 2), ValidationError);
  assert.throws(() => M.lineEvidence(1, 1, 0), ValidationError);
});
