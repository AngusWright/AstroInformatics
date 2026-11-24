"use strict";
const { test } = require("node:test");
const assert = require("node:assert/strict");
const G = require("./gibbs.js");

function near(actual, expected, tolerance, message = "") {
  assert.ok(Math.abs(actual - expected) <= tolerance, `${message}: ${actual} != ${expected} (+/- ${tolerance})`);
}

function integrate(f, lo, hi, n = 12000) {
  const step = (hi - lo) / n;
  let sum = f(lo) + f(hi);
  for (let i = 1; i < n; i++) sum += (i % 2 ? 4 : 2) * f(lo + i * step);
  return sum * step / 3;
}

test("replay, prefix invariance and sequential full-conditional draws", () => {
  for (const family of ["gaussian", "student", "mixture", "poisson"]) {
    const chain = G.sampleChain({ family, sweeps: 25 });
    assert.deepEqual(chain.states, G.sampleChain({ family, sweeps: 25 }).states);
    assert.deepEqual(chain.states.slice(0, 21), G.sampleChain({ family, sweeps: 10 }).states);
    const rng = G.createRng(chain.seed);
    for (let update = 1; update < chain.states.length; update++) {
      const view = G.frame(chain, update), j = (update - 1) % 2;
      assert.equal(view.coordinate, j);
      assert.equal(view.after[1 - j], view.before[1 - j]);
      assert.equal(view.after[j], chain.target.conditional(view.before, j).draw(rng));
      if (family === "poisson") assert.ok(view.after.every(x => x > 0));
    }
    assert.notDeepEqual(chain.states, G.sampleChain({ family, seed: 1, sweeps: 25 }).states);
  }
});

test("incomplete sweeps and burn-in never enter summaries", () => {
  const chain = G.sampleChain({ sweeps: 5 });
  assert.equal(G.frame(chain, 0).hasDraw, false);
  assert.deepEqual(G.samples(chain, 1), []);
  assert.deepEqual(G.samples(chain, 3), [chain.states[2]]);
  assert.deepEqual(G.samples(chain, 7, 1), [chain.states[4], chain.states[6]]);
  assert.deepEqual(G.samples(chain, 10, 5), []);
  near(G.summarize([1, 2, 3]).sd, 1, 1e-14);
  assert.deepEqual(G.acf([1, 2, 3], 2), [1, 0, -0.5]);
  assert.throws(() => G.acf([1, 1, 1]), G.ValidationError);
});

test("conditional densities agree with normalised slices of every joint", () => {
  for (const family of ["gaussian", "student", "mixture", "poisson"]) {
    const target = G.createTarget({ family, rho: -0.7, weight: 0.3, separation: 5 });
    for (const j of [0, 1]) {
      const fixed = family === "poisson" ? [2, 1.3] : [1.2, -0.8];
      const d = target.conditional(fixed, j);
      const lo = target.positive ? 0 : -120, hi = target.positive ? 100 : 120;
      near(integrate(d.density, lo, hi), 1, 1e-6, `${family} density integral`);
      const xs = family === "poisson" ? [0.8, 2, 4] : [-2, 0, 2];
      const joint = x => j === 0 ? target.logDensity(x, fixed[1]) : target.logDensity(fixed[0], x);
      for (const x of xs) {
        near(Math.log(d.density(x)) - Math.log(d.density(xs[0])),
          joint(x) - joint(xs[0]), 1e-10, `${family} conditional/joint identity`);
      }
    }
  }
});

test("Gaussian and Student-t conditional parameters and marginal moments", () => {
  const settings = { location1: 1, location2: -2, scale1: 2, scale2: 3, rho: 0.6, df: 5 };
  const normal = G.createTarget(settings).conditional([4, 1], 0);
  near(normal.mean, 2.2, 1e-12);
  near(normal.sd, 1.6, 1e-12);
  const target = G.createTarget({ ...settings, family: "student" });
  const conditional = target.conditional([4, 1], 0);
  near(conditional.mean, normal.mean, 1e-12);
  near(conditional.sd, 1.6 * Math.sqrt(6 / 4), 1e-12);
  near(target.sd[0], 2 * Math.sqrt(5 / 3), 1e-12);
  for (const family of ["gaussian", "student"]) {
    const chain = G.sampleChain({ ...settings, family, sweeps: 60000, seed: 401 });
    const kept = G.samples(chain, chain.states.length - 1, 1000);
    for (const j of [0, 1]) {
      const summary = G.summarize(kept.map(s => s[j]));
      near(summary.mean, chain.target.mean[j], 0.1, `${family} stationary mean`);
      near(summary.sd, chain.target.sd[j], 0.13, `${family} stationary SD`);
    }
  }
});

test("mixture weights use the held coordinate; zero separation collapses", () => {
  const target = G.createTarget({ family: "mixture", rho: 0, weight: 0.3, separation: 4 });
  const q = 1 / (1 + Math.exp(-(Math.log(0.3 / 0.7) - 4)));
  near(target.conditional([0, 1], 0).mean, q * -2 + (1 - q) * 2, 1e-12);
  assert.equal(target.conditional([100, 1], 0).mean, target.conditional([-100, 1], 0).mean);
  assert.ok(target.conditional([0, -1000], 0).mean < -1.99);
  assert.ok(target.conditional([0, 1000], 0).mean > 1.99);
  const normal = G.createTarget({ rho: -0.98 }), collapsed = G.createTarget({
    family: "mixture", separation: 0, weight: 0.2, rho: -0.98
  });
  for (const j of [0, 1]) for (const x of [-3, 0, 2]) {
    near(collapsed.conditional([2, -3], j).density(x), normal.conditional([2, -3], j).density(x), 1e-12);
    near(collapsed.marginal(x, j), normal.marginal(x, j), 1e-12);
  }
  near(collapsed.logDensity(0, 1), normal.logDensity(0, 1), 1e-12);
  const chain = G.sampleChain({ family: "mixture", separation: 2, weight: 0.3, rho: 0, sweeps: 60000, seed: 23 });
  const kept = G.samples(chain, 120000, 1000);
  for (const j of [0, 1]) {
    const summary = G.summarize(kept.map(s => s[j]));
    near(summary.mean, chain.target.mean[j], 0.05, "mixture stationary mean");
    near(summary.sd, chain.target.sd[j], 0.04, "mixture stationary SD");
  }
});

test("Gamma rate parameterisation and Poisson posterior stationary moments", () => {
  const rng = G.createRng(902);
  const draws = Array.from({ length: 60000 }, () => rng.gamma(0.7, 2));
  const summary = G.summarize(draws);
  near(summary.mean, 0.7 / 2, 0.008);
  near(summary.sd, Math.sqrt(0.7) / 2, 0.01);
  const target = G.createTarget({ family: "poisson" });
  near(target.conditional([3, 2], 0).mean, 32 / 12, 1e-12);
  near(target.conditional([3, 2], 1).mean, 4 / 4, 1e-12);
  const chain = G.sampleChain({ family: "poisson", sweeps: 60000, seed: 27 });
  const kept = G.samples(chain, 120000, 1000);
  // Integrate beta out: p(lambda|x) is proportional to lambda^31 exp(-10 lambda)/(1+lambda)^4.
  const marginal = x => x === 0 ? 0 : Math.exp(31 * Math.log(x) - 10 * x - 4 * Math.log1p(x));
  const mass = integrate(marginal, 0, 25);
  const mean = integrate(x => x * marginal(x), 0, 25) / mass;
  near(G.summarize(kept.map(s => s[0])).mean, mean, 0.015);
  const betaMean = integrate(x => 4 / (1 + x) * marginal(x), 0, 25) / mass;
  near(G.summarize(kept.map(s => s[1])).mean, betaMean, 0.015);
});

test("invalid settings fail explicitly", () => {
  for (const settings of [
    { family: "unknown" }, { rho: 1 }, { scale1: 0 }, { location2: NaN },
    { family: "student", df: 2 }, { family: "mixture", weight: 1 },
    { family: "mixture", separation: -1 }, { family: "poisson", n: 0 },
    { family: "poisson", total: -1 }, { family: "poisson", hyper_rate: 0 },
    { family: "poisson", start_beta: 0 }, { start1: Infinity }, { seed: -1 }, { sweeps: 1.5 }
  ]) assert.throws(() => G.sampleChain(settings), G.ValidationError);
  const chain = G.sampleChain({ sweeps: 2 });
  assert.throws(() => G.frame(chain, 5), G.ValidationError);
  assert.throws(() => G.samples(chain, 2, -1), G.ValidationError);
});
