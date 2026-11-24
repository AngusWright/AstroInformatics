(function (root, factory) {
  "use strict";
  const model = factory(typeof module === "object" && module.exports ? require("../TeachingTools/statistics.js") : root.Astro);
  if (typeof module === "object" && module.exports) module.exports = model;
  else root.CentralLimitModel = model;
})(globalThis, function (A) {
  "use strict";
  const populations = {
    exponential: { mean: 1, variance: 1, limits: [-0.2, 6], label: "Exponential observations" },
    mixture: { mean: 0.25, variance: 2.1875, limits: [-4, 9], label: "95% N(0,1) + 5% N(5,1)" },
    cauchy: { mean: null, variance: null, limits: [-10, 10], label: "Standard Cauchy observations" }
  };
  function validate(p) {
    A.check(Object.hasOwn(populations, p.population), "Choose a supported observation population.");
    A.check(Number.isInteger(p.seed) && p.seed >= 0 && p.seed <= 2147483647, "Seed must be an integer from 0 to 2147483647.");
    A.check(Number.isInteger(p.n) && p.n >= 1 && p.n <= 100, "Measurements per average must be an integer from 1 to 100.");
    A.check(Number.isInteger(p.repetitions) && p.repetitions >= 200 && p.repetitions <= 3000,
      "Use between 200 and 3000 independent averaging experiments.");
  }
  function populationDensity(population, x) {
    A.check(Object.hasOwn(populations, population), "Choose a supported observation population.");
    if (population === "exponential") return x < 0 ? 0 : Math.exp(-x);
    if (population === "mixture") return 0.95 * A.normalPdf(x) + 0.05 * A.normalPdf(x, 5, 1);
    return 1 / (Math.PI * (1 + x * x));
  }
  function prepare(seed, population) {
    validate({ seed, population, n: 1, repetitions: 200 });
    const rng = A.createRng(seed);
    const draw = population === "exponential" ? () => -Math.log1p(-rng.uniform()) :
      population === "mixture" ? () => (rng.uniform() < 0.05 ? 5 : 0) + rng.normal() :
        () => Math.tan(Math.PI * (rng.uniform() - 0.5));
    const observations = Array.from({ length: 8000 }, draw);
    const sums = new Float64Array(3000 * 100);
    for (let experiment = 0; experiment < 3000; experiment++) {
      let sum = 0;
      for (let j = 0; j < 100; j++) {
        sum += draw();
        sums[experiment * 100 + j] = sum;
      }
    }
    return { seed, population, observations, sums };
  }
  function compute(parameters, prepared) {
    validate(parameters);
    const data = prepared || prepare(parameters.seed, parameters.population);
    A.check(data.seed === parameters.seed && data.population === parameters.population && data.sums.length === 300000,
      "Cached observations do not match these settings.");
    return { data, parameters: { ...parameters }, population: populations[parameters.population] };
  }
  function atFrame(result, frame = 100) {
    A.check(Number.isInteger(frame) && frame >= 0 && frame <= 100, "Animation frame must be from 0 to 100.");
    const p = result.parameters, population = result.population;
    const n = Math.max(1, Math.round(1 + frame / 100 * (p.n - 1)));
    const means = Array.from({ length: p.repetitions }, (_, i) => result.data.sums[i * 100 + n - 1] / n);
    const finite = population.variance !== null;
    const statistic = finite ? means.map(x => Math.sqrt(n) * (x - population.mean) / Math.sqrt(population.variance)) : means;
    return { n, means, statistic, finite, meanSd: finite ? Math.sqrt(population.variance / n) : null,
      empiricalIqr: A.quantile(means, 0.75) - A.quantile(means, 0.25) };
  }
  return { populations, validate, populationDensity, prepare, compute, atFrame };
});
