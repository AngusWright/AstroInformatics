(function (root) {
  "use strict";
  const A = typeof module !== "undefined" && module.exports ? require("../TeachingTools/statistics.js") : root.Astro;
  const DEFAULTS = Object.freeze({
    model: "normal", no_data: false, log_n: 1, sample_mean: 1, sigma: 2,
    prior_mean: -1, prior_sd: 1, success_fraction: 0.7, beta_a: 2, beta_b: 5,
    count_mean: 4, gamma_mean: 2, shape: 3, gamma_a: 3, gamma_b: 2
  });
  const MODELS = Object.freeze({
    normal: {
      description: "Infer the mean of Normal observations with known variance. Prior and posterior are Normal.",
      rule: "theta ~ Normal(b, variance=s0^2)\nd = 1/s0^2; tau = 1/sigma^2\n\ntheta | x ~ Normal(mean=(d*b + n*tau*xbar)/(d+n*tau),\n                   variance=1/(d+n*tau))",
      weight: "Prior weight = d/(d + n/sigma^2), where d = 1/s0^2."
    },
    binomial: {
      description: "Infer the probability of success from s successes in n independent trials. Prior and posterior are Beta.",
      rule: "theta ~ Beta(alpha, beta)\ns | theta ~ Binomial(n, theta)\n\ntheta | s ~ Beta(alpha+s, beta+n-s)",
      weight: "Prior weight = (alpha+beta)/(alpha+beta+n). Prior concentration acts as an effective sample size for this mean update."
    },
    poisson: {
      description: "Infer a Poisson rate from equal unit-exposure counts. Prior and posterior are Gamma.",
      rule: "theta ~ Gamma(shape=alpha, rate=beta)\nx_i | theta ~ Poisson(theta)\n\ntheta | x ~ Gamma(shape=alpha+sum(x), rate=beta+n)",
      weight: "Prior weight = beta/(beta+n). The prior rate beta acts as prior exposure in units of one observation."
    },
    gamma: {
      description: "Infer the rate of Gamma observations with known shape k. Prior and posterior are Gamma.",
      rule: "theta ~ Gamma(shape=alpha, rate=beta)\nx_i | theta ~ Gamma(shape=k, rate=theta)\n\ntheta | x ~ Gamma(shape=alpha+n*k, rate=beta+sum(x))",
      weight: "Prior weight = beta/(beta+sum(x)). The data estimate of rate is n*k/sum(x), not the sample mean."
    }
  });
  function roundEven(x) {
    const lo = Math.floor(x);
    return x - lo === 0.5 ? (lo % 2 === 0 ? lo : lo + 1) : Math.round(x);
  }
  function update(model, n, statistic, prior1, prior2, { sigma = 1, shape = 1 } = {}) {
    A.check(Object.hasOwn(MODELS, model), "Choose a supported conjugate pair.");
    A.check(Number.isSafeInteger(n) && n >= 0, "Sample size must be a non-negative integer.");
    A.check(Number.isFinite(statistic) && Number.isFinite(prior1) && Number.isFinite(prior2) && prior2 > 0,
      "Use finite sufficient statistics and prior parameters; prior spread or rate must be positive.");
    if (model !== "normal") {
      A.check(prior1 > 0 && statistic >= 0, "Prior shape must be positive; sufficient statistic must be non-negative.");
      A.check(n !== 0 || statistic === 0, "With no observations the count or sum must be zero.");
    }
    if (model === "binomial") A.check(Number.isSafeInteger(statistic) && statistic <= n, "Successes must be an integer between zero and n.");
    if (model === "poisson") A.check(Number.isSafeInteger(statistic), "Count total must be an integer.");
    if (model === "gamma") {
      A.check(Number.isFinite(shape) && shape > 0 && (n === 0 || statistic > 0),
        "Gamma observation shape and nonempty data sum must be positive.");
    }
    let prior, posterior, likelihood = null, estimate = null, weight, description, parameters;
    if (model === "normal") {
      A.check(Number.isFinite(sigma) && sigma > 0, "Known observation SD must be positive.");
      prior = A.normal(prior1, prior2);
      const d = 1 / prior2 ** 2, precision = d + n / sigma ** 2;
      weight = d / precision;
      parameters = [weight * prior1 + (1 - weight) * statistic, 1 / Math.sqrt(precision)];
      posterior = n === 0 ? prior : A.normal(...parameters);
      if (n > 0) { likelihood = A.normal(statistic, sigma / Math.sqrt(n)); estimate = statistic; }
      description = `n = ${n}; observed sample mean = ${statistic}; known sigma = ${sigma}.`;
    } else if (model === "binomial") {
      prior = A.beta(prior1, prior2);
      parameters = [prior1 + statistic, prior2 + n - statistic];
      posterior = n === 0 ? prior : A.beta(...parameters);
      weight = (prior1 + prior2) / (prior1 + prior2 + n);
      if (n > 0) { likelihood = A.beta(statistic + 1, n - statistic + 1); estimate = statistic / n; }
      description = `n = ${n} trials; ${statistic} successes; ${n - statistic} failures.`;
    } else {
      prior = A.gamma(prior1, prior2);
      parameters = model === "poisson" ? [prior1 + statistic, prior2 + n] : [prior1 + n * shape, prior2 + statistic];
      posterior = n === 0 ? prior : A.gamma(...parameters);
      weight = prior2 / parameters[1];
      if (n > 0) {
        likelihood = model === "poisson" ? A.gamma(statistic + 1, n) : A.gamma(n * shape + 1, statistic);
        estimate = model === "poisson" ? statistic / n : n * shape / statistic;
      }
      description = model === "poisson" ? `n = ${n} unit-exposure observations; total count = ${statistic}.` :
        `n = ${n}; sum(x) = ${Number(statistic.toPrecision(6))}; known observation shape k = ${shape}.`;
    }
    return { model, n, statistic, prior, posterior, likelihood, estimate, prior_weight: weight,
      posterior_parameters: parameters, data_description: n === 0 ? "n = 0: no observed data. Posterior equals prior." : description,
      weight_description: MODELS[model].weight };
  }
  function analyse(settings = {}, nOverride) {
    const p = { ...DEFAULTS, ...settings };
    A.check(Object.hasOwn(MODELS, p.model), "Choose a supported conjugate pair.");
    A.check(Number.isFinite(p.log_n) && p.log_n >= 0 && p.log_n <= 4, "Choose a sample-size exponent between zero and four.");
    const n = nOverride === undefined ? (p.no_data ? 0 : roundEven(10 ** p.log_n)) : nOverride;
    let statistic, prior;
    if (p.model === "normal") { statistic = p.sample_mean; prior = [p.prior_mean, p.prior_sd]; }
    if (p.model === "binomial") {
      A.check(Number.isFinite(p.success_fraction) && p.success_fraction >= 0 && p.success_fraction <= 1, "Success fraction must be between zero and one.");
      statistic = roundEven(n * p.success_fraction); prior = [p.beta_a, p.beta_b];
    }
    if (p.model === "poisson") {
      A.check(Number.isFinite(p.count_mean) && p.count_mean >= 0, "Observed counts per observation must be non-negative.");
      statistic = roundEven(n * p.count_mean); prior = [p.gamma_a, p.gamma_b];
    }
    if (p.model === "gamma") {
      A.check(Number.isFinite(p.gamma_mean) && p.gamma_mean > 0, "Gamma observation mean must be positive.");
      statistic = n * p.gamma_mean; prior = [p.gamma_a, p.gamma_b];
    }
    return update(p.model, n, statistic, ...prior, { sigma: p.sigma, shape: p.shape });
  }
  function roleCurve(settings = {}) {
    const ns = [...new Set(Array.from({ length: 121 }, (_, i) => roundEven(10 ** (4 * i / 120))))];
    return ns.map(n => {
      const r = analyse(settings, n);
      return { n, mean: r.posterior.mean, lower: r.posterior.quantile(0.025), upper: r.posterior.quantile(0.975),
        prior_weight: r.prior_weight, estimate: r.estimate };
    });
  }
  const api = { DEFAULTS, MODELS, update, analyse, roleCurve };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  root.ConjugateModel = api;
})(globalThis);
