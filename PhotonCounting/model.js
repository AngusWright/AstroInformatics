(function (root, factory) {
  "use strict";
  const api = factory(typeof module === "object" && module.exports ? require("../TeachingTools/statistics.js") : root.Astro);
  if (typeof module === "object" && module.exports) module.exports = api;
  root.PhotonCountingModel = api;
})(globalThis, function (A) {
  "use strict";
  const defaults = { seed: 42, rate: 2, exposure: 5, gaussian: "show" };
  const STEPS = 60, REPEATS = 3000, MAX_LAMBDA = 160;
  function validate(p) {
    A.check(Number.isInteger(p.seed) && p.seed >= 0 && p.seed <= 2147483647, "Enter an integer seed between 0 and 2147483647.");
    A.check(Number.isFinite(p.rate) && p.rate >= 0 && p.rate <= 8, "Photon rate must be between 0 and 8 per second.");
    A.check(Number.isFinite(p.exposure) && p.exposure >= 0.5 && p.exposure <= 20, "Exposure must be between 0.5 and 20 seconds.");
    A.check(["show", "hide"].includes(p.gaussian), "Choose whether to show the Gaussian approximation.");
  }
  function checkLambda(lambda) {
    A.check(Number.isFinite(lambda) && lambda >= 0 && lambda <= MAX_LAMBDA, "This lesson supports Poisson means between 0 and 160.");
  }
  function poissonPmf(k, lambda) {
    checkLambda(lambda);
    A.check(Number.isInteger(k), "Poisson counts must be integers.");
    if (k < 0) return 0;
    return lambda === 0 ? Number(k === 0) : Math.exp(-lambda + k * Math.log(lambda) - A.logGamma(k + 1));
  }
  function poissonTable(lambda) {
    checkLambda(lambda);
    const upper = Math.max(5, Math.ceil(lambda + 10 * Math.sqrt(lambda) + 10));
    const pmf = [Math.exp(-lambda)], cdf = [pmf[0]];
    let sum = pmf[0], compensation = 0;
    for (let k = 1; k <= upper; k++) {
      pmf.push(pmf[k - 1] * lambda / k);
      const add = pmf[k] - compensation, next = sum + add;
      compensation = (next - sum) - add; sum = next;
      cdf.push(Math.min(1, sum));
    }
    return { pmf, cdf };
  }
  function fromTable(probability, table) {
    let lo = 0, hi = table.cdf.length - 1;
    A.check(table.cdf[hi] >= probability, "The requested Poisson tail exceeds numerical precision.");
    while (lo < hi) {
      const mid = Math.floor((lo + hi) / 2);
      if (table.cdf[mid] >= probability) hi = mid;
      else lo = mid + 1;
    }
    return lo;
  }
  function poissonQuantile(probability, lambda) {
    checkLambda(lambda);
    A.check(Number.isFinite(probability) && probability >= 0 && probability <= 1, "Probability must be between zero and one.");
    if (lambda === 0 || probability === 0) return 0;
    if (probability === 1) return Infinity;
    return fromTable(probability, poissonTable(lambda));
  }
  function gaussianMass(k, lambda) {
    checkLambda(lambda);
    A.check(Number.isInteger(k), "Gaussian comparison bins must be centred on integers.");
    if (lambda === 0) return null;
    const sd = Math.sqrt(lambda);
    return A.normalCdf(k + 0.5, lambda, sd) - A.normalCdf(k - 0.5, lambda, sd);
  }
  function prepare(seed) {
    const rng = A.createRng(seed);
    const uniforms = Array.from({ length: REPEATS }, () => rng.uniform());
    const unitArrivals = [];
    let time = -Math.log(rng.uniform());
    while (time <= MAX_LAMBDA) { unitArrivals.push(time); time -= Math.log(rng.uniform()); }
    return { seed, uniforms, unitArrivals };
  }
  function compute(parameters, prepared) {
    const p = { ...defaults, ...parameters };
    validate(p);
    const base = prepared || prepare(p.seed);
    A.check(base.seed === p.seed && base.uniforms.length === REPEATS, "The cached photon experiment must match the seed.");
    const allArrivals = p.rate === 0 ? [] : base.unitArrivals.map(time => time / p.rate).filter(time => time <= p.exposure);
    const frames = Array.from({ length: STEPS + 1 }, (_, frame) => {
      const exposure = p.exposure * frame / STEPS, lambda = p.rate * exposure;
      const distribution = poissonTable(lambda);
      const counts = base.uniforms.map(u => fromTable(u, distribution));
      const summary = A.summarize(counts);
      const lower = Math.min(-3, Math.floor(lambda - 4 * Math.sqrt(lambda)));
      const upper = Math.max(5, Math.max(...counts), Math.ceil(lambda + 4 * Math.sqrt(lambda)));
      const integers = Array.from({ length: upper - lower + 1 }, (_, i) => lower + i);
      const empirical = Array(upper + 1).fill(0);
      counts.forEach(count => { empirical[count] += 1 / REPEATS; });
      return { exposure, lambda, arrivals: allArrivals.filter(time => time <= exposure), counts,
        mean: summary.mean, variance: summary.sd ** 2, sd: Math.sqrt(lambda),
        relativeSd: lambda > 0 ? 1 / Math.sqrt(lambda) : null,
        negativeMass: lambda > 0 ? A.normalCdf(-0.5, lambda, Math.sqrt(lambda)) : null,
        lower, upper, integers, empirical,
        exact: integers.map(k => k < 0 ? 0 : distribution.pmf[k] || 0),
        gaussian: lambda > 0 ? integers.map(k => gaussianMass(k, lambda)) : null };
    });
    return { base, frames, allArrivals, rate: p.rate, exposure: p.exposure, steps: STEPS, repeats: REPEATS };
  }
  function atFrame(result, frame) {
    A.check(Number.isInteger(frame) && frame >= 0 && frame <= STEPS, "Choose an integer photon frame from 0 to 60.");
    return result.frames[frame];
  }
  return { defaults, STEPS, REPEATS, MAX_LAMBDA, validate, poissonPmf, poissonTable, poissonQuantile, gaussianMass,
    prepare, compute, atFrame };
});
