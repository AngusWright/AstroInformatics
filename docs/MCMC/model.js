/* Numerical model and sampler. No DOM access or automatic sampling. */
(function (root, factory) {
  const api = factory(typeof module === 'object' && module.exports
    ? require('../TeachingTools/statistics.js') : root.Astro);
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.MCMCModel = api;
})(typeof globalThis === 'object' ? globalThis : this, function (A) {
  'use strict';
  const LOG2PI = Math.log(2 * Math.PI);
  const check = A.check;
  const finite = x => typeof x === 'number' && Number.isFinite(x);
  const integer = (x, lo, hi) => Number.isInteger(x) && x >= lo && x <= hi;
  const normalLog = (x, mean, sd) => -0.5 * (LOG2PI + ((x - mean) / sd) ** 2) - Math.log(sd);
  const logAdd = (a, b) => {
    if (a === -Infinity) return b;
    if (b === -Infinity) return a;
    const m = Math.max(a, b);
    return m + Math.log(Math.exp(a - m) + Math.exp(b - m));
  };

  function defaults(family = 'gaussian') {
    return {
      family, theta: family === 'gaussian' ? [0, 0] : family === 'binomial' ? [0.3, 5] : [7, 2],
      observations: 10, dataSd: 1.2, rho: 0.8,
      bounds: family === 'gaussian' ? [[-3.5, 3.5], [-3.5, 3.5]]
        : family === 'binomial' ? [[0, 1], [1, 10]] : [[0, 10], [0, 10]],
      prior: 'uniform', proposal: 'randomwalk', proposalSd: 0.2, globalWeight: 0.2,
      iterations: 1000, chains: 4, seed: 666, starts: null
    };
  }

  function validate(c) {
    check(c && ['gaussian', 'binomial', 'beta', 'gamma'].includes(c.family), 'Choose a supported model.');
    check(Array.isArray(c.theta) && c.theta.length === 2 && c.theta.every(finite), 'Both data-generating parameters must be finite.');
    check(integer(c.observations, 1, 10000), 'Number of observations must be an integer from 1 to 10,000.');
    check(Array.isArray(c.bounds) && c.bounds.length === 2 && c.bounds.every(b =>
      Array.isArray(b) && b.length === 2 && b.every(finite) && b[0] < b[1] && b[1] - b[0] <= 200),
    'Each prior range needs finite, increasing endpoints (width at most 200).');
    check(['uniform', 'gaussian'].includes(c.prior), 'Choose a uniform or bounded Gaussian prior.');
    check(['randomwalk', 'independence', 'mixture'].includes(c.proposal), 'Choose a supported proposal.');
    check(finite(c.proposalSd) && c.proposalSd > 0 && c.proposalSd <= 100, 'Proposal SD must be positive and at most 100.');
    check(finite(c.globalWeight) && c.globalWeight > 0 && c.globalWeight < 1, 'The global mixture weight must lie strictly between 0 and 1.');
    check(integer(c.iterations, 1, 10000), 'Iterations must be an integer from 1 to 10,000.');
    check(integer(c.chains, 1, 8), 'Run between 1 and 8 chains.');
    check(integer(c.seed, 0, 2147483647), 'Seed must be an integer from 0 to 2,147,483,647.');
    if (c.family === 'gaussian') {
      check(finite(c.dataSd) && c.dataSd > 0, 'Data SD must be positive.');
      check(finite(c.rho) && Math.abs(c.rho) < 1, 'Correlation must lie strictly between −1 and 1; endpoints give a singular covariance.');
    } else if (c.family === 'binomial') {
      check(c.theta[0] >= 0 && c.theta[0] <= 1 && integer(c.theta[1], 1, 100), 'Binomial data need probability in [0, 1] and integer draws from 1 to 100.');
      check(c.bounds[0][0] >= 0 && c.bounds[0][1] <= 1, 'The probability prior must be contained in [0, 1].');
      check(c.bounds[1].every(x => integer(x, 1, 100)), 'The draw-count prior needs integer endpoints from 1 to 100.');
    } else {
      check(c.theta.every(x => x > 0), 'Beta and Gamma data-generating shapes/rates must be strictly positive.');
      check(c.bounds.every(b => b[0] >= 0), 'Shape and rate prior ranges cannot extend below zero.');
    }
    return c;
  }

  function simulateData(config) {
    const c = validate(config), rng = A.createRng(c.seed), values = [];
    for (let i = 0; i < c.observations; i++) {
      let value;
      if (c.family === 'gaussian') {
        const z = rng.normal(), w = rng.normal();
        value = [c.theta[0] + c.dataSd * z,
          c.theta[1] + c.dataSd * (c.rho * z + Math.sqrt(1 - c.rho * c.rho) * w)];
      } else if (c.family === 'binomial') {
        value = 0;
        for (let j = 0; j < c.theta[1]; j++) if (rng.uniform() < c.theta[0]) value++;
      } else if (c.family === 'beta') {
        const x = rng.gamma(c.theta[0], 1), y = rng.gamma(c.theta[1], 1);
        value = x / (x + y);
        check(finite(value) && value > 0 && value < 1, 'The simulated Beta data reached floating-point limits. Use less extreme shapes.');
      } else {
        value = rng.gamma(c.theta[0], c.theta[1]);
        check(finite(value) && value > 0, 'The simulated Gamma data reached floating-point limits. Use less extreme parameters.');
      }
      values.push(value);
    }
    return values;
  }

  function createTarget(config, suppliedData) {
    const c = validate(config);
    const data = suppliedData === undefined ? simulateData(c) : suppliedData.map(x => Array.isArray(x) ? x.slice() : x);
    check(data.length === c.observations, 'The supplied data count must equal the observation count.');
    const n = data.length, binomial = c.family === 'binomial';
    const midpoint = c.bounds.map(b => (b[0] + b[1]) / 2);
    const priorSd = c.bounds.map(b => (b[1] - b[0]) / 4);
    const discreteValues = binomial
      ? Array.from({ length: c.bounds[1][1] - c.bounds[1][0] + 1 }, (_, i) => c.bounds[1][0] + i) : null;
    const gaussianNorm = c.bounds.map((b, d) => binomial && d === 1
      ? Math.log(discreteValues.reduce((s, k) => s + Math.exp(normalLog(k, midpoint[d], priorSd[d])), 0))
      : Math.log(A.normalCdf(b[1], midpoint[d], priorSd[d]) - A.normalCdf(b[0], midpoint[d], priorSd[d])));
    let mean, scatter, total = 0, sumLog = 0, sumLog1m = 0, maximum = 0;
    const counts = new Map();
    if (c.family === 'gaussian') {
      check(data.every(v => Array.isArray(v) && v.length === 2 && v.every(finite)), 'Gaussian observations must be finite coordinate pairs.');
      mean = [0, 1].map(d => data.reduce((s, v) => s + v[d], 0) / n);
      scatter = data.reduce((s, v) => {
        const x = v[0] - mean[0], y = v[1] - mean[1];
        return s + x * x - 2 * c.rho * x * y + y * y;
      }, 0);
    } else {
      for (const x of data) {
        check(finite(x) && (binomial ? Number.isInteger(x) && x >= 0
          : c.family === 'beta' ? x > 0 && x < 1 : x > 0), 'An observation is outside the model support.');
        total += x;
        if (binomial) {
          maximum = Math.max(maximum, x);
          counts.set(x, (counts.get(x) || 0) + 1);
        } else {
          sumLog += Math.log(x);
          if (c.family === 'beta') sumLog1m += Math.log1p(-x);
        }
      }
      if (binomial) check(maximum <= c.bounds[1][1], 'The prior permits no draw count large enough for the observed data. Increase its upper limit.');
    }
    function inSupport(x) {
      return Array.isArray(x) && x.length === 2 && x.every((v, d) =>
        finite(v) && v >= c.bounds[d][0] && v <= c.bounds[d][1]) &&
        (binomial ? Number.isInteger(x[1]) && x[1] >= maximum
          : c.family === 'gaussian' || x.every(v => v > 0));
    }
    function logPrior(x) {
      if (!inSupport(x)) return -Infinity;
      if (c.prior === 'uniform') return -Math.log(c.bounds[0][1] - c.bounds[0][0])
        - Math.log(binomial ? discreteValues.length : c.bounds[1][1] - c.bounds[1][0]);
      return x.reduce((s, v, d) => s + normalLog(v, midpoint[d], priorSd[d]) - gaussianNorm[d], 0);
    }
    function logLikelihood(x) {
      if (!inSupport(x)) return -Infinity;
      const [a, b] = x;
      if (c.family === 'gaussian') {
        const dx = a - mean[0], dy = b - mean[1], variance = c.dataSd ** 2;
        return -n * (LOG2PI + Math.log(variance) + 0.5 * Math.log1p(-(c.rho ** 2)))
          - (scatter + n * (dx * dx - 2 * c.rho * dx * dy + dy * dy)) / (2 * variance * (1 - c.rho ** 2));
      }
      if (binomial) {
        if ((a === 0 && total > 0) || (a === 1 && total < n * b)) return -Infinity;
        let comb = 0;
        for (const [k, count] of counts) comb += count * (A.logGamma(b + 1) - A.logGamma(k + 1) - A.logGamma(b - k + 1));
        return comb + (total === 0 ? 0 : total * Math.log(a))
          + (n * b === total ? 0 : (n * b - total) * Math.log1p(-a));
      }
      if (c.family === 'beta') return n * (A.logGamma(a + b) - A.logGamma(a) - A.logGamma(b))
        + (a - 1) * sumLog + (b - 1) * sumLog1m;
      return n * (a * Math.log(b) - A.logGamma(a)) + (a - 1) * sumLog - b * total;
    }
    function logDensity(x) {
      const prior = logPrior(x);
      return prior === -Infinity ? -Infinity : prior + logLikelihood(x);
    }
    return { config: c, data, binomial, maximum, midpoint, priorSd, discreteValues,
      logPrior, logLikelihood, logDensity, inSupport,
      labels: c.family === 'gaussian' ? ['Mean x', 'Mean y']
        : binomial ? ['Probability p', 'Draw count N'] : c.family === 'beta' ? ['Shape α', 'Shape β'] : ['Shape k', 'Rate r'],
      dataSummary: c.family === 'gaussian' ? `Observed mean: (${mean.map(x => x.toPrecision(4)).join(', ')}).`
        : binomial ? `Observed mean: ${(total / n).toPrecision(4)}; largest count: ${maximum}.`
          : `Observed mean: ${(total / n).toPrecision(4)}.`
    };
  }

  function localLogQ(target, from, to) {
    const sd = target.config.proposalSd;
    if (target.binomial) {
      const jump = to[1] - from[1];
      return Number.isInteger(jump) && Math.abs(jump) <= 1
        ? normalLog(to[0], from[0], sd) - Math.log(3) : -Infinity;
    }
    return normalLog(to[0], from[0], sd) + normalLog(to[1], from[1], sd);
  }
  function globalLogQ(target, to) {
    return to.reduce((s, v, d) => {
      if (target.binomial && d === 1) return Number.isInteger(v) &&
        v >= target.config.bounds[1][0] && v <= target.config.bounds[1][1]
        ? s - Math.log(target.discreteValues.length) : -Infinity;
      return s + normalLog(v, target.midpoint[d], (target.config.bounds[d][1] - target.config.bounds[d][0]) / 3);
    }, 0);
  }
  function logProposal(target, from, to) {
    const c = target.config;
    if (c.proposal === 'randomwalk') return localLogQ(target, from, to);
    if (c.proposal === 'independence') return globalLogQ(target, to);
    return logAdd(Math.log1p(-c.globalWeight) + localLogQ(target, from, to),
      Math.log(c.globalWeight) + globalLogQ(target, to));
  }
  function propose(target, from, rng) {
    const c = target.config;
    const global = c.proposal === 'independence' || (c.proposal === 'mixture' && rng.uniform() < c.globalWeight);
    const proposal = global
      ? c.bounds.map((b, d) => target.binomial && d === 1
        ? b[0] + Math.floor(rng.uniform() * target.discreteValues.length)
        : target.midpoint[d] + (b[1] - b[0]) / 3 * rng.normal())
      : [from[0] + c.proposalSd * rng.normal(),
        from[1] + (target.binomial ? Math.floor(3 * rng.uniform()) - 1 : c.proposalSd * rng.normal())];
    return { proposal, component: global ? 'Global independence' : 'Local random walk' };
  }
  function transition(target, from, rng) {
    check(Number.isFinite(target.logDensity(from)), 'The current state must have positive, finite target density.');
    const { proposal, component } = propose(target, from, rng);
    const logTargetBefore = target.logDensity(from), logTargetProposal = target.logDensity(proposal);
    const logQForward = logProposal(target, from, proposal), logQReverse = logProposal(target, proposal, from);
    const logRatio = logTargetProposal === -Infinity ? -Infinity
      : logTargetProposal - logTargetBefore + logQReverse - logQForward;
    check(!Number.isNaN(logRatio), 'The log acceptance ratio is not numerically defined. Use less extreme parameters.');
    const logAlpha = Math.min(0, logRatio), logU = Math.log(rng.uniform()), accepted = logU < logAlpha;
    return { before: from.slice(), proposal, state: accepted ? proposal.slice() : from.slice(),
      accepted, component, logTargetBefore, logTargetProposal, logQForward, logQReverse, logRatio, logAlpha, logU };
  }
  function run(config, suppliedData) {
    const c = validate(config), target = createTarget(c, suppliedData);
    const rng = A.createRng((c.seed ^ 0x1e3779b9) & 0x7fffffff);
    if (c.starts !== null && c.starts !== undefined) {
      check(Array.isArray(c.starts) && c.starts.length === c.chains, 'Supply one starting pair per chain.');
      check(c.starts.every(s => Number.isFinite(target.logDensity(s))), 'Every starting state must have positive target density.');
    }
    const chains = [];
    for (let j = 0; j < c.chains; j++) {
      const start = c.starts ? c.starts[j].slice() : c.bounds.map((b, d) => {
        if (target.binomial && d === 1) {
          const low = Math.max(b[0], target.maximum);
          return low + Math.floor(rng.uniform() * (b[1] - low + 1));
        }
        return b[0] + rng.uniform() * (b[1] - b[0]);
      });
      check(Number.isFinite(target.logDensity(start)), 'The starting state has non-finite target density. Check the prior and data.');
      const states = [start], moves = [];
      for (let t = 0; t < c.iterations; t++) {
        const move = transition(target, states[t], rng);
        moves.push(move);
        states.push(move.state);
      }
      chains.push({ states, moves });
    }
    return { config: c, target, data: target.data, chains };
  }

  function retained(result, cursor, burn, thin) {
    check(integer(cursor, 0, result.config.iterations), 'Playback position is outside the cached chain.');
    check(integer(burn, 0, result.config.iterations), 'Burn-in must be an integer from zero to the iteration count.');
    check(integer(thin, 1, 10000), 'Thinning must be a positive integer no greater than 10,000.');
    return result.chains.map(chain => {
      const values = [];
      for (let t = burn + 1; t <= cursor; t += thin) values.push(chain.states[t]);
      return values;
    });
  }

  function maximize(fn, lo, hi) {
    const ratio = (Math.sqrt(5) - 1) / 2;
    let a = lo, b = hi, x = b - ratio * (b - a), y = a + ratio * (b - a);
    let fx = fn(x), fy = fn(y);
    for (let i = 0; i < 55; i++) {
      if (fx > fy) { b = y; y = x; fy = fx; x = b - ratio * (b - a); fx = fn(x); }
      else { a = x; x = y; fx = fy; y = a + ratio * (b - a); fy = fn(y); }
    }
    const points = [lo, hi, x, y];
    return points.reduce((best, v) => fn(v) > fn(best) ? v : best, x);
  }
  function modeOf(target) {
    const bounds = target.config.bounds;
    if (target.binomial) {
      let best = null, bestLog = -Infinity;
      for (const k of target.discreteValues) {
        if (k < target.maximum) continue;
        const p = maximize(x => target.logDensity([x, k]), ...bounds[0]);
        const value = target.logDensity([p, k]);
        if (value > bestLog) { best = [p, k]; bestLog = value; }
      }
      return best;
    }
    let state = target.midpoint.slice();
    for (let i = 0; i < 100; i++) {
      const before = state.slice();
      for (let d = 0; d < 2; d++) state[d] = maximize(v => {
        const x = state.slice(); x[d] = v; return target.logDensity(x);
      }, ...bounds[d]);
      if (Math.hypot(...state.map((v, d) => v - before[d])) < 1e-9) break;
    }
    return state;
  }
  function axisGrid(lo, hi, center, scale, cells) {
    const edges = [lo];
    for (const [a, b, side] of [[lo, center, -1], [center, hi, 1]]) {
      if (b - a < (hi - lo) * 1e-12) continue;
      const distance = b - a, bend = Math.asinh(distance / scale);
      for (let i = 1; i <= cells; i++) {
        const t = side === -1 ? 1 - i / cells : i / cells;
        edges.push(i === cells ? b : center + side * scale * Math.sinh(t * bend));
      }
    }
    edges[edges.length - 1] = hi;
    return { edges, centers: edges.slice(1).map((v, i) => (v + edges[i]) / 2),
      widths: edges.slice(1).map((v, i) => v - edges[i]) };
  }

  function reference(target, cells = 85) {
    check(integer(cells, 15, 250), 'Reference grid resolution must be an integer from 15 to 250.');
    const bounds = target.config.bounds, mode = modeOf(target);
    check(mode && Number.isFinite(target.logDensity(mode)), 'No finite posterior density was found inside the prior.');
    const scales = bounds.map(b => (b[1] - b[0]) / 8);
    const probe = mode.map((v, d) => Math.min(bounds[d][1] - (bounds[d][1] - bounds[d][0]) * 1e-5,
      Math.max(bounds[d][0] + (bounds[d][1] - bounds[d][0]) * 1e-5, v)));
    if (!target.binomial) {
      const h = bounds.map(b => (b[1] - b[0]) * 1e-6);
      const f = (a, b) => target.logDensity([probe[0] + a * h[0], probe[1] + b * h[1]]);
      const f0 = f(0, 0);
      const h00 = -(f(1, 0) - 2 * f0 + f(-1, 0)) / (h[0] * h[0]);
      const h11 = -(f(0, 1) - 2 * f0 + f(0, -1)) / (h[1] * h[1]);
      const h01 = -(f(1, 1) - f(1, -1) - f(-1, 1) + f(-1, -1)) / (4 * h[0] * h[1]);
      const determinant = h00 * h11 - h01 * h01;
      if (h00 > 0 && h11 > 0 && determinant > 0) {
        scales[0] = Math.sqrt(h11 / determinant); scales[1] = Math.sqrt(h00 / determinant);
      }
    } else {
      const h = (bounds[0][1] - bounds[0][0]) * 1e-6, x = probe[0], k = mode[1];
      const curvature = -(target.logDensity([x + h, k]) - 2 * target.logDensity([x, k])
        + target.logDensity([x - h, k])) / (h * h);
      if (curvature > 0) scales[0] = 1 / Math.sqrt(curvature);
    }
    for (let d = 0; d < 2; d++) {
      const width = bounds[d][1] - bounds[d][0];
      if (!(target.binomial && d === 1) &&
          Math.min(mode[d] - bounds[d][0], bounds[d][1] - mode[d]) < width * 1e-6) {
        const low = mode.slice(), high = mode.slice(), h = width * 1e-6;
        low[d] = probe[d] - h; high[d] = probe[d] + h;
        const slope = (target.logDensity(high) - target.logDensity(low)) / (2 * h);
        // A boundary peak can have an exponential width much smaller than its curvature scale.
        if (Number.isFinite(slope) && slope !== 0) scales[d] = Math.min(scales[d], 1 / Math.abs(slope));
      }
      scales[d] = Math.max(width * 1e-12, Math.min(width, scales[d]));
    }
    const axes = bounds.map((b, d) => target.binomial && d === 1
      ? { centers: target.discreteValues, edges: [...target.discreteValues.map(k => k - 0.5), b[1] + 0.5],
        widths: target.discreteValues.map(() => 1) }
      : axisGrid(b[0], b[1], mode[d], scales[d], cells));
    const nx = axes[0].centers.length, ny = axes[1].centers.length;
    const logs = [], marginal = [new Array(nx).fill(0), new Array(ny).fill(0)];
    let maxLog = -Infinity;
    for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) {
      const value = target.logDensity([axes[0].centers[i], axes[1].centers[j]]);
      logs.push(value); maxLog = Math.max(maxLog, value);
    }
    let normalizer = 0;
    const mass = logs.map((v, k) => {
      const i = k % nx, j = Math.floor(k / nx);
      const w = Math.exp(v - maxLog) * axes[0].widths[i] * axes[1].widths[j];
      normalizer += w; marginal[0][i] += w; marginal[1][j] += w;
      return w;
    });
    check(normalizer > 0 && Number.isFinite(normalizer), 'The posterior reference could not be normalised.');
    const moments = marginal.map((weights, d) => {
      let mean = 0;
      weights.forEach((v, i) => { mean += v / normalizer * axes[d].centers[i]; });
      let variance = 0;
      weights.forEach((v, i) => { variance += v / normalizer * (axes[d].centers[i] - mean) ** 2; });
      return { mean, sd: Math.sqrt(variance) };
    });
    return { axes, mode, nx, ny, logs, maxLog, mass: mass.map(v => v / normalizer), moments,
      marginal: marginal.map((weights, d) => weights.map((v, i) => v / normalizer / axes[d].widths[i])) };
  }

  function diagnostics(chains, dimension) {
    const values = chains.map(chain => chain.map(x => x[dimension]));
    const pooled = values.flat(), summary = pooled.length >= 2 ? A.summarize(pooled)
      : { mean: pooled.length ? pooled[0] : null, sd: null };
    const length = Math.min(...values.map(x => x.length)), maxLag = Math.min(40, length - 1);
    let acf = [];
    if (length >= 3 && values.every(v => A.summarize(v).sd > 0)) {
      const individual = values.map(v => A.acf(v, maxLag));
      acf = Array.from({ length: maxLag + 1 }, (_, lag) =>
        individual.reduce((s, a) => s + a[lag], 0) / individual.length);
    }
    let ess = null;
    if (acf.length >= 3 && acf.every(finite)) {
      let sum = 0;
      for (let k = 1; k + 1 < acf.length; k += 2) {
        const pair = acf[k] + acf[k + 1];
        if (pair <= 0) break;
        sum += pair;
      }
      ess = pooled.length / Math.max(1, 1 + 2 * sum);
    }
    let rhat = null;
    const n = Math.floor(length / 2);
    if (values.length >= 2 && n >= 2) {
      const halves = values.flatMap(v => [v.slice(0, n), v.slice(-n)]).map(A.summarize);
      const means = halves.map(v => v.mean), W = halves.reduce((s, v) => s + v.sd ** 2, 0) / halves.length;
      const B = n * A.summarize(means).sd ** 2;
      if (W > 0) rhat = Math.sqrt(((n - 1) * W / n + B / n) / W);
    }
    return { ...summary, acf, ess, rhat, count: pooled.length };
  }
  return { defaults, validate, simulateData, createTarget, logProposal, propose, transition, run, retained,
    reference, diagnostics };
});
