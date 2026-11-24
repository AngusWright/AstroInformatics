"use strict";
const { test } = require("node:test");
const assert = require("node:assert/strict");
const A = require("../TeachingTools/statistics.js");
const C = require("./model.js");
function near(actual, expected, tolerance = 1e-10) {
  assert.ok(Math.abs(actual - expected) <= tolerance, `${actual} != ${expected}`);
}

test("all four exact conjugate updates and mean decompositions", () => {
  const normal = C.update("normal", 4, 3, 1, 2, { sigma: 2 });
  near(normal.posterior.mean, 2.6); near(normal.posterior.sd ** 2, 0.8);
  near(normal.prior_weight, 0.2);
  const binomial = C.update("binomial", 10, 7, 2, 5);
  near(binomial.posterior.mean, 9 / 17);
  near(binomial.posterior.quantile(0.025), 0.29877689908317751);
  const poisson = C.update("poisson", 10, 40, 3, 2);
  near(poisson.posterior.mean, 43 / 12); near(poisson.posterior.sd ** 2, 43 / 144);
  const gamma = C.update("gamma", 5, 10, 3, 2, { shape: 3 });
  near(gamma.posterior.mean, 18 / 12); near(gamma.posterior.sd ** 2, 18 / 144);
  for (const r of [normal, binomial, poisson, gamma]) {
    near(r.posterior.mean, r.prior_weight * r.prior.mean + (1 - r.prior_weight) * r.estimate);
  }
});

test("posterior agrees with the prior times actual observation likelihood", () => {
  const data = [1, 2, 3, 6];
  const cases = [
    ["normal", 3, theta => data.reduce((s, x) => s + Math.log(A.normalPdf(x, theta, 2)), 0)],
    ["binomial", 3, theta => 3 * Math.log(theta) + Math.log1p(-theta)],
    ["poisson", 6, theta => 6 * Math.log(theta) - 4 * theta],
    ["gamma", 6, theta => 12 * Math.log(theta) - 6 * theta]
  ];
  for (const [model, statistic, logLikelihood] of cases) {
    const r = C.update(model, 4, statistic, 2, 3, { sigma: 2, shape: 3 });
    const xs = model === "binomial" ? [0.2, 0.4, 0.7] : [0.5, 1, 2];
    const products = xs.map(x => Math.log(r.prior.density(x)) + logLikelihood(x));
    const post = xs.map(x => Math.log(r.posterior.density(x)));
    for (let i = 1; i < xs.length; i++) near(products[i] - products[0], post[i] - post[0]);
  }
});

test("sequential batches equal one update using the combined sufficient statistic", () => {
  for (const [model, firstN, firstStat, secondN, secondStat] of [
    ["normal", 2, 1, 2, 5], ["binomial", 3, 2, 5, 4],
    ["poisson", 3, 6, 5, 10], ["gamma", 3, 6, 5, 10]
  ]) {
    const options = { sigma: 2, shape: 3 };
    const first = C.update(model, firstN, firstStat, 2, 3, options);
    const second = C.update(model, secondN, secondStat, ...first.posterior_parameters, options);
    const combinedStat = model === "normal" ? (firstN * firstStat + secondN * secondStat) / (firstN + secondN) : firstStat + secondStat;
    const combined = C.update(model, firstN + secondN, combinedStat, 2, 3, options);
    near(second.posterior.mean, combined.posterior.mean);
    near(second.posterior.sd, combined.posterior.sd);
  }
});

test("no data, boundary counts, fixed observations and asymptotic influence", () => {
  for (const model of Object.keys(C.MODELS)) {
    const empty = C.analyse({ model, no_data: true });
    assert.equal(empty.posterior, empty.prior);
    assert.equal(empty.likelihood, null); assert.equal(empty.estimate, null);
    near(empty.prior_weight, 1);
    const original = C.analyse({ model });
    const changed = C.analyse({ model, prior_mean: 4, beta_a: 12, gamma_a: 10 });
    assert.equal(changed.statistic, original.statistic);
    for (const x of [0.2, 0.5, 0.8]) near(changed.likelihood.density(x), original.likelihood.density(x));
  }
  assert.equal(C.analyse({ model: "binomial", success_fraction: 0.25 }).statistic, 2);
  for (const fraction of [0, 1]) {
    const r = C.analyse({ model: "binomial", log_n: 4, success_fraction: fraction, beta_a: 0.5, beta_b: 0.5 });
    assert.ok(r.posterior.mean > 0 && r.posterior.mean < 1);
    assert.ok([0.001, 0.999].every(p => Number.isFinite(r.posterior.quantile(p))));
  }
  const zero = C.analyse({ model: "poisson", log_n: 4, count_mean: 0, gamma_a: 0.5, gamma_b: 0.1 });
  assert.ok(zero.posterior.mean > 0);
  const weak = C.analyse({ prior_sd: 1e6 });
  near(weak.posterior.mean, 1); near(weak.posterior.sd ** 2, 0.4);
  const curve = C.roleCurve();
  assert.ok(curve.at(-1).prior_weight < curve[0].prior_weight);
  assert.ok(curve.every((r, i) => i === 0 || r.upper - r.lower < curve[i - 1].upper - curve[i - 1].lower));
});

test("invalid sufficient statistics and priors are rejected explicitly", () => {
  for (const args of [
    ["normal", 1, 2, 0, 0], ["normal", 1, 2, 0, 1, { sigma: 0 }],
    ["binomial", 10, 11, 1, 1], ["binomial", 10, 2.5, 1, 1],
    ["poisson", 0, 1, 1, 1], ["poisson", 2, 0.5, 1, 1],
    ["gamma", 10, 0, 1, 1], ["gamma", 10, 1, 1, 1, { shape: 0 }]
  ]) assert.throws(() => C.update(...args), A.ValidationError);
  assert.throws(() => C.analyse({ model: "unknown" }), A.ValidationError);
});
