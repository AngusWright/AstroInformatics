(function (root) {
  "use strict";

  const DEFAULTS = Object.freeze({
    family: "gaussian", location1: 0, location2: 0, scale1: 1, scale2: 1,
    rho: 0.8, df: 5, separation: 4, weight: 0.5,
    n: 10, total: 30, a: 2, hyper_shape: 2, hyper_rate: 1,
    sweeps: 100, seed: 666, start1: -2.5, start2: 2.5, start_lambda: 1, start_beta: 1
  });
  const LOG_2PI = Math.log(2 * Math.PI);
  const { ValidationError, check: requireValue, createRng, logGamma, normal, gamma, summarize, acf } =
    typeof module !== "undefined" && module.exports ? require("../TeachingTools/statistics.js") : root.Astro;

  function student(location, scale, df) {
    const sd = scale * Math.sqrt(df / (df - 2));
    const constant = logGamma((df + 1) / 2) - logGamma(df / 2) -
      Math.log(df * Math.PI) / 2 - Math.log(scale);
    return {
      label: `Student-t: df ${df}, location ${location.toPrecision(4)}, scale ${scale.toPrecision(4)}`,
      mean: location, sd, range: [location - 6 * sd, location + 6 * sd],
      focus: [{ location, scale }],
      draw: rng => location + scale * rng.normal() / Math.sqrt(rng.gamma(df / 2, 0.5) / df),
      density: x => Math.exp(constant - (df + 1) / 2 * Math.log1p(((x - location) / scale) ** 2 / df))
    };
  }

  function mixture(locations, scale, weight) {
    const parts = locations.map(location => normal(location, scale));
    const mean = weight * locations[0] + (1 - weight) * locations[1];
    return {
      label: `Conditional mixture: weights ${weight.toPrecision(3)} / ${(1 - weight).toPrecision(3)}; ` +
        `means ${locations.map(x => x.toPrecision(4)).join(" / ")}; component SD ${scale.toPrecision(4)}`,
      mean, sd: Math.sqrt(scale ** 2 + weight * (1 - weight) * (locations[1] - locations[0]) ** 2),
      range: [Math.min(...locations) - 4.5 * scale, Math.max(...locations) + 4.5 * scale],
      focus: parts.flatMap(part => part.focus),
      draw: rng => parts[rng.uniform() < weight ? 0 : 1].draw(rng),
      density: x => weight * parts[0].density(x) + (1 - weight) * parts[1].density(x)
    };
  }

  function logistic(x) {
    return x >= 0 ? 1 / (1 + Math.exp(-x)) : Math.exp(x) / (1 + Math.exp(x));
  }

  function createTarget(settings = {}) {
    const p = { ...DEFAULTS, ...settings };
    requireValue(["gaussian", "student", "mixture", "poisson"].includes(p.family), "Unknown target family.");
    if (p.family === "poisson") {
      requireValue(Number.isInteger(p.n) && p.n >= 1 && Number.isInteger(p.total) && p.total >= 0,
        "Use a positive integer observation count and a non-negative integer total.");
      requireValue([p.a, p.hyper_shape, p.hyper_rate].every(Number.isFinite) &&
        p.a >= 1 && p.hyper_shape >= 1 && p.hyper_rate > 0, "Gamma shapes must be at least 1; rate must be positive.");
      return {
        family: p.family, label: "Poisson with Gamma-Gamma hierarchy", axes: ["lambda", "beta"],
        positive: true, mean: null, sd: null, limits: null, marginal: null,
        conditional: (state, j) => j === 0 ? gamma(p.a + p.total, p.n + state[1]) :
          gamma(p.a + p.hyper_shape, p.hyper_rate + state[0]),
        logDensity: (x, y) => x <= 0 || y <= 0 ? -Infinity :
          (p.total + p.a - 1) * Math.log(x) + (p.a + p.hyper_shape - 1) * Math.log(y) -
          p.n * x - x * y - p.hyper_rate * y
      };
    }
    const location = [p.location1, p.location2];
    const scale = [p.scale1, p.scale2];
    requireValue(location.every(Number.isFinite) && scale.every(x => Number.isFinite(x) && x > 0) &&
      Number.isFinite(p.rho) && Math.abs(p.rho) < 1, "Locations must be finite, scales positive, and |rho| < 1.");
    if (p.family === "mixture") {
      requireValue(Number.isFinite(p.separation) && p.separation >= 0 &&
        Number.isFinite(p.weight) && p.weight > 0 && p.weight < 1,
      "Separation must be non-negative; population weight must be strictly between 0 and 1.");
      const parts = [-1, 1].map(sign => createTarget({
        ...p, family: "gaussian",
        location1: p.location1 + sign * p.separation * p.scale1 / 2,
        location2: p.location2 + sign * p.separation * p.scale2 / 2
      }));
      return {
        family: p.family, label: "Two-population Gaussian mixture", axes: ["theta1", "theta2"], positive: false,
        mean: location.map((x, j) => x + (1 - 2 * p.weight) * p.separation * scale[j] / 2),
        sd: scale.map(s => s * Math.sqrt(1 + p.weight * (1 - p.weight) * p.separation ** 2)),
        limits: [0, 1].map(j => [Math.min(...parts.map(part => part.limits[j][0])),
          Math.max(...parts.map(part => part.limits[j][1]))]),
        conditional(state, j) {
          const other = 1 - j;
          const q = logistic(Math.log(p.weight / (1 - p.weight)) -
            p.separation * ((state[other] - location[other]) / scale[other]));
          return mixture(parts.map(part => part.conditional(state, j).mean),
            scale[j] * Math.sqrt(1 - p.rho ** 2), q);
        },
        marginal: (x, j) => p.weight * parts[0].marginal(x, j) + (1 - p.weight) * parts[1].marginal(x, j),
        logDensity(x, y) {
          const a = Math.log(p.weight) + parts[0].logDensity(x, y);
          const b = Math.log1p(-p.weight) + parts[1].logDensity(x, y);
          const peak = Math.max(a, b);
          return peak === -Infinity ? peak : peak + Math.log(Math.exp(a - peak) + Math.exp(b - peak));
        }
      };
    }
    requireValue(p.family !== "student" || (Number.isFinite(p.df) && p.df > 2),
      "Student-t degrees of freedom must exceed 2.");
    const marginals = location.map((x, j) => p.family === "student" ?
      student(x, scale[j], p.df) : normal(x, scale[j]));
    return {
      family: p.family, label: p.family === "student" ? "Bivariate Student-t posterior" : "Bivariate Gaussian posterior",
      axes: ["theta1", "theta2"], positive: false, mean: location, sd: marginals.map(d => d.sd),
      limits: marginals.map(d => d.range), marginal: (x, j) => marginals[j].density(x),
      conditional(state, j) {
        const other = 1 - j;
        const residual = (state[other] - location[other]) / scale[other];
        const centre = location[j] + p.rho * scale[j] * residual;
        const width = scale[j] * Math.sqrt(1 - p.rho ** 2);
        return p.family === "student" ?
          student(centre, width * Math.sqrt((p.df + residual ** 2) / (p.df + 1)), p.df + 1) : normal(centre, width);
      },
      logDensity(x, y) {
        const u = (x - location[0]) / scale[0];
        const v = (y - location[1]) / scale[1];
        const q = (u * u - 2 * p.rho * u * v + v * v) / (1 - p.rho ** 2);
        const determinant = Math.log(scale[0]) + Math.log(scale[1]) + Math.log1p(-(p.rho ** 2)) / 2;
        return p.family === "student" ?
          logGamma((p.df + 2) / 2) - logGamma(p.df / 2) - Math.log(p.df * Math.PI) -
            determinant - (p.df + 2) / 2 * Math.log1p(q / p.df) :
          -LOG_2PI - determinant - q / 2;
      }
    };
  }

  function sampleChain(settings = {}) {
    const p = { ...DEFAULTS, ...settings };
    const target = createTarget(p);
    requireValue(Number.isInteger(p.sweeps) && p.sweeps >= 1, "Sweep count must be a positive integer.");
    const rng = createRng(p.seed);
    const start = target.positive ? [p.start_lambda, p.start_beta] : [p.start1, p.start2];
    requireValue(start.every(Number.isFinite) && (!target.positive || start.every(x => x > 0)),
      "Enter finite starting coordinates; Poisson starts must both be positive.");
    const states = [start];
    for (let update = 0; update < 2 * p.sweeps; update++) {
      const j = update % 2;
      const next = states[update].slice();
      next[j] = target.conditional(next, j).draw(rng);
      requireValue(Number.isFinite(next[j]) && (!target.positive || next[j] > 0),
        "A draw exceeded numerical precision. Try less extreme settings.");
      states.push(next);
    }
    return { target, states, sweeps: p.sweeps, seed: p.seed };
  }

  function frame(chain, update) {
    requireValue(Number.isInteger(update) && update >= 0 && update < chain.states.length, "Update is outside this chain.");
    const coordinate = update === 0 ? 0 : (update - 1) % 2;
    const before = chain.states[Math.max(0, update - 1)];
    return { update, coordinate, before, after: chain.states[update], hasDraw: update > 0,
      distribution: chain.target.conditional(before, coordinate) };
  }

  function samples(chain, update, burn = 0) {
    frame(chain, update);
    requireValue(Number.isInteger(burn) && burn >= 0, "Burn-in must be a non-negative integer.");
    const result = [];
    for (let sweep = burn + 1; sweep <= Math.floor(update / 2); sweep++) result.push(chain.states[2 * sweep]);
    return result;
  }

  const api = { DEFAULTS, ValidationError, createRng, logGamma, normal, student, gamma, mixture,
    createTarget, sampleChain, frame, samples, summarize, acf };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  root.Gibbs = api;
})(globalThis);
