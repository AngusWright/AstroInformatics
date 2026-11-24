'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const A = require('../TeachingTools/statistics.js');
const M = require('./model.js');
const near = (actual, expected, tolerance = 1e-10) =>
  assert.ok(Math.abs(actual - expected) <= tolerance, `${actual} should be within ${tolerance} of ${expected}`);
function configuration(family, overrides = {}) { return { ...M.defaults(family), ...overrides }; }
function fakeRng(normals, uniforms) {
  let ni = 0, ui = 0;
  return {
    normal() { assert.ok(ni < normals.length, 'Unexpected normal draw'); return normals[ni++]; },
    uniform() { assert.ok(ui < uniforms.length, 'Unexpected uniform draw'); return uniforms[ui++]; }
  };
}

test('all four likelihoods retain their parameter-dependent normalisers', () => {
  let c = configuration('gaussian', { observations: 2, dataSd: 1.3, rho: 0.6 });
  let t = M.createTarget(c, [[0.2, 0.8], [-0.4, 0.1]]);
  const state = [0.3, -0.2];
  const expected = t.data.reduce((s, z) => {
    const dx = z[0] - state[0], dy = z[1] - state[1];
    return s - Math.log(2 * Math.PI * 1.3 ** 2 * Math.sqrt(1 - 0.6 ** 2))
      - (dx * dx - 1.2 * dx * dy + dy * dy) / (2 * 1.3 ** 2 * (1 - 0.6 ** 2));
  }, 0);
  near(t.logLikelihood(state), expected);

  c = configuration('binomial', { observations: 3 });
  t = M.createTarget(c, [0, 2, 1]);
  near(t.logLikelihood([0.3, 3]), Math.log(0.7 ** 3) + Math.log(3 * 0.3 ** 2 * 0.7) + Math.log(3 * 0.3 * 0.7 ** 2));
  assert.equal(t.logDensity([0.3, 2.5]), -Infinity);
  assert.equal(t.logDensity([0.3, 1]), -Infinity);
  assert.equal(t.logDensity([0, 3]), -Infinity);
  assert.equal(t.logDensity([1, 3]), -Infinity);

  t = M.createTarget(configuration('beta', { observations: 2 }), [0.2, 0.7]);
  near(t.logLikelihood([2, 3]), Math.log(12 * 0.2 * 0.8 ** 2) + Math.log(12 * 0.7 * 0.3 ** 2));
  assert.equal(t.logDensity([0, 3]), -Infinity);
  t = M.createTarget(configuration('gamma', { observations: 2 }), [0.4, 2]);
  near(t.logLikelihood([2, 3]), Math.log(9 * 0.4) - 1.2 + Math.log(9 * 2) - 6);
  assert.equal(t.logDensity([2, 0]), -Infinity);
});

test('a rejection records the old state and invalid boundary proposals are not retried', () => {
  const c = configuration('gaussian', { observations: 1, dataSd: 1, rho: 0, proposalSd: 1, bounds: [[-2, 2], [-2, 2]] });
  const t = M.createTarget(c, [[0, 0]]);
  const rejected = M.transition(t, [0, 0], fakeRng([1, 0], [0.9]));
  near(rejected.logRatio, -0.5);
  assert.equal(rejected.accepted, false);
  assert.deepEqual(rejected.proposal, [1, 0]);
  assert.deepEqual(rejected.state, [0, 0]);
  const accepted = M.transition(t, [0, 0], fakeRng([1, 0], [0.1]));
  assert.equal(accepted.accepted, true);
  assert.deepEqual(accepted.state, [1, 0]);
  const outside = M.transition(t, [0, 0], fakeRng([3, 0], [1e-100]));
  assert.equal(outside.logAlpha, -Infinity);
  assert.equal(outside.accepted, false);
  assert.deepEqual(outside.state, [0, 0]);

  const b = M.createTarget(configuration('binomial', { observations: 1 }), [2]);
  const invalidN = M.transition(b, [0.4, 2], fakeRng([0], [0.1, 0.01]));
  assert.deepEqual(invalidN.proposal, [0.4, 1]);
  assert.equal(invalidN.accepted, false);
  assert.deepEqual(invalidN.state, [0.4, 2]);
  const zero = M.createTarget(configuration('binomial', { observations: 1 }), [0]);
  assert.ok(Number.isFinite(zero.logDensity([0, 1])));
  assert.equal(zero.logDensity([-0.01, 1]), -Infinity);
});

test('asymmetric and mixture proposals use full forward and reverse densities', () => {
  const c = configuration('gaussian', { observations: 1, dataSd: 1, rho: 0,
    bounds: [[-2, 2], [-2, 2]], proposal: 'mixture', globalWeight: 0.3, proposalSd: 0.4 });
  const t = M.createTarget(c, [[0, 0]]), from = [0, 0], to = [1.2, -0.4];
  const local = (x, y) => A.normalPdf(y[0], x[0], 0.4) * A.normalPdf(y[1], x[1], 0.4);
  const global = y => A.normalPdf(y[0], 0, 4 / 3) * A.normalPdf(y[1], 0, 4 / 3);
  const forward = 0.7 * local(from, to) + 0.3 * global(to);
  const reverse = 0.7 * local(to, from) + 0.3 * global(from);
  near(M.logProposal(t, from, to), Math.log(forward));
  near(M.logProposal(t, to, from), Math.log(reverse));
  const move = M.transition(t, from, fakeRng([0.9, -0.3], [0.1, 0.5]));
  near(move.logRatio, t.logDensity(to) - t.logDensity(from) + Math.log(reverse / forward));
  assert.notEqual(move.logQForward, move.logQReverse);
  assert.notEqual(move.logRatio, t.logDensity(to) - t.logDensity(from));

  const independent = M.createTarget({ ...c, proposal: 'independence' }, [[0, 0]]);
  near(M.logProposal(independent, from, to), Math.log(global(to)));
  near(M.logProposal(independent, to, from), Math.log(global(from)));
});

test('invalid settings fail explicitly, including singular covariance and impossible data', () => {
  for (const change of [{ rho: 1 }, { rho: -1 }, { dataSd: 0 }, { observations: 1.5 },
    { proposalSd: 0 }, { bounds: [[1, -1], [0, 2]] }, { seed: -1 }, { globalWeight: 1 }]) {
    assert.throws(() => M.run(configuration('gaussian', change)), A.ValidationError);
  }
  assert.throws(() => M.run(configuration('beta', { theta: [0, 2] })), A.ValidationError);
  assert.throws(() => M.run(configuration('gamma', { theta: [2, -1] })), A.ValidationError);
  assert.throws(() => M.run(configuration('binomial', { theta: [0.3, 2.5] })), A.ValidationError);
  assert.throws(() => M.createTarget(configuration('binomial', { observations: 1 }), [11]), A.ValidationError);
  assert.throws(() => M.run(configuration('gaussian', { chains: 1, starts: [[99, 99]] })), A.ValidationError);
});

test('seeded complete runs replay exactly and sampler changes leave observations unchanged', () => {
  const c = configuration('gamma', { iterations: 60, chains: 2 });
  const first = M.run(c), second = M.run(c);
  assert.deepEqual(first.data, second.data);
  assert.deepEqual(first.chains, second.chains);
  const changedProposal = M.run({ ...c, proposal: 'mixture', proposalSd: 0.8 });
  assert.deepEqual(first.data, changedProposal.data);
  assert.notDeepEqual(first.chains, changedProposal.chains);
  assert.notDeepEqual(first.data, M.run({ ...c, seed: c.seed + 1 }).data);
  const before = JSON.stringify(first.chains);
  const kept = M.retained(first, 17, 4, 3);
  assert.deepEqual(kept[0], [5, 8, 11, 14, 17].map(t => first.chains[0].states[t]));
  assert.deepEqual(M.retained(first, 3, 4, 3), [[], []]);
  assert.equal(JSON.stringify(first.chains), before);
  assert.throws(() => M.retained(first, 61, 0, 1), A.ValidationError);
});

test('full-support quadrature resolves narrow correlated Gaussian posteriors', () => {
  const c = configuration('gaussian', { observations: 10000, dataSd: 1.2, rho: 0.98 });
  const target = M.createTarget(c, Array.from({ length: c.observations }, () => [0.12, -0.3]));
  const ref = M.reference(target);
  near(ref.mass.reduce((s, x) => s + x, 0), 1, 1e-12);
  for (let d = 0; d < 2; d++) {
    near(ref.moments[d].mean, [0.12, -0.3][d], 0.0001);
    near(ref.moments[d].sd, 0.012, 0.00008);
  }
  near(ref.axes[0].edges[0], c.bounds[0][0]);
  near(ref.axes[0].edges.at(-1), c.bounds[0][1]);
});

test('reference resolves a narrow boundary posterior rather than fitting the visited region', () => {
  const c = configuration('gaussian', { observations: 1000, dataSd: 1, rho: 0, bounds: [[0, 1], [0, 1]] });
  const target = M.createTarget(c, Array.from({ length: c.observations }, () => [2, 2]));
  const ref = M.reference(target);
  const sd = 1 / Math.sqrt(1000), z = -1 / sd, ratio = A.normalPdf(z) / A.normalCdf(z);
  const mean = 2 - sd * ratio, targetSd = sd * Math.sqrt(1 - z * ratio - ratio * ratio);
  for (const moment of ref.moments) {
    near(moment.mean, mean, 0.000005);
    near(moment.sd, targetSd, 0.000005);
  }
});

test('bounded Gaussian priors are applied to continuous and integer parameters', () => {
  const c = configuration('gaussian', { observations: 1, prior: 'gaussian', rho: 0, dataSd: 0.5,
    bounds: [[-5, 5], [-5, 5]] });
  const target = M.createTarget(c, [[0.4, -0.2]]), ref = M.reference(target);
  const variance = 1 / (1 / 0.5 ** 2 + 1 / 2.5 ** 2);
  for (let d = 0; d < 2; d++) {
    near(ref.moments[d].mean, [0.4, -0.2][d] * variance / 0.5 ** 2, 0.0001);
    near(ref.moments[d].sd, Math.sqrt(variance), 0.0001);
  }
  const b = M.createTarget(configuration('binomial', { observations: 1, prior: 'gaussian' }), [0]);
  const center = 0.5, step = 0.001;
  const numerical = Array.from({ length: 1000 }, (_, i) => (i + 0.5) * step)
    .reduce((sum, p) => sum + b.discreteValues.reduce((s, k) => s + Math.exp(b.logPrior([p, k])) * step, 0), 0);
  near(numerical, 1, 0.000001);
  assert.ok(b.logPrior([center, 5]) > b.logPrior([center, 1]));
});

test('all three fixed Gaussian proposal kernels sample the correct stationary target', () => {
  for (const proposal of ['randomwalk', 'independence', 'mixture']) {
    const c = configuration('gaussian', {
      observations: 1, dataSd: 0.7, rho: 0.4, bounds: [[-5, 5], [-5, 5]],
      proposal, proposalSd: 0.8, globalWeight: 0.35, iterations: 10000, chains: 4, seed: 127
    });
    const result = M.run(c, [[0.4, -0.2]]), kept = M.retained(result, 10000, 1500, 1);
    for (let d = 0; d < 2; d++) {
      const s = M.diagnostics(kept, d);
      near(s.mean, [0.4, -0.2][d], 0.055);
      near(s.sd, 0.7, 0.045);
    }
    assert.ok(result.chains.every(chain => chain.moves.some(m => !m.accepted) && chain.moves.some(m => m.accepted)));
  }
});

test('mixed continuous/discrete Binomial target agrees with exact integrated probabilities', () => {
  const c = configuration('binomial', { observations: 3, bounds: [[0, 1], [1, 5]],
    iterations: 10000, chains: 4, proposal: 'mixture', proposalSd: 0.12, globalWeight: 0.35, seed: 125 });
  const result = M.run(c, [0, 1, 2]), ref = M.reference(result.target);
  const weights = [2, 3, 4, 5].map(k => Math.exp(2 * Math.log(k) + Math.log(k - 1) - Math.log(2)
    + A.logGamma(4) + A.logGamma(3 * k - 2) - A.logGamma(3 * k + 2)));
  const sum = weights.reduce((s, x) => s + x, 0);
  const nMean = weights.reduce((s, x, i) => s + (i + 2) * x / sum, 0);
  const pMean = weights.reduce((s, x, i) => s + 4 / (3 * (i + 2) + 2) * x / sum, 0);
  near(ref.moments[0].mean, pMean, 0.0001);
  near(ref.moments[1].mean, nMean, 0.0001);
  const kept = M.retained(result, c.iterations, 1000, 1);
  near(M.diagnostics(kept, 0).mean, pMean, 0.015);
  near(M.diagnostics(kept, 1).mean, nMean, 0.08);
  assert.ok(result.chains.every(chain => chain.states.every(x => Number.isInteger(x[1]) && x[1] >= 2)));
});

test('Beta and Gamma samplers agree with full bounded posterior references', () => {
  for (const [family, data] of [['beta', [0.2, 0.45, 0.6]], ['gamma', [1, 2, 4]]]) {
    const c = configuration(family, { observations: 3, bounds: [[0, 5], [0, 5]], iterations: 10000,
      chains: 4, proposal: 'mixture', proposalSd: 0.7, globalWeight: 0.35, seed: 12 });
    const result = M.run(c, data), ref = M.reference(result.target);
    const kept = M.retained(result, c.iterations, 1000, 1);
    for (let d = 0; d < 2; d++) {
      near(M.diagnostics(kept, d).mean, ref.moments[d].mean, 0.08);
      near(M.diagnostics(kept, d).sd, ref.moments[d].sd, 0.08);
    }
    assert.ok(result.chains.every(chain => chain.states.every(x => x.every(v => v > 0 && v <= 5))));
  }
});

test('diagnostics do not concatenate chains and handle absent or constant samples', () => {
  const chains = [[0, 1, 0, 1], [10, 11, 10, 11]].map(c => c.map(x => [x, x]));
  const diagnostics = M.diagnostics(chains, 0);
  near(diagnostics.acf[1], -0.75);
  assert.ok(diagnostics.rhat > 5);
  assert.deepEqual(M.diagnostics([[], []], 0), { mean: null, sd: null, acf: [], ess: null, rhat: null, count: 0 });
  assert.equal(M.diagnostics([[[2, 2]]], 0).mean, 2);
  assert.equal(M.diagnostics([[[2, 2]]], 0).sd, null);
  const constant = M.diagnostics([Array.from({ length: 10 }, () => [2, 2])], 0);
  assert.deepEqual(constant.acf, []);
  assert.equal(constant.ess, null);
});
