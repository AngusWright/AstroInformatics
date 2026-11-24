/* Numerical model only: no DOM access. */
(function (root, factory) {
  "use strict";
  const api = factory(typeof module === "object" && module.exports
    ? require("../TeachingTools/statistics.js") : root.Astro);
  if (typeof module === "object" && module.exports) module.exports = api;
  else root.ShiftModel = api;
})(typeof globalThis !== "undefined" ? globalThis : this, function (Astro) {
  "use strict";
  const { check, createRng } = Astro;
  const Z95 = 1.959963984540054;
  const grid = (a, b, n) => Array.from({ length: n }, (_, i) => a + (b - a) * i / (n - 1));
  const sum = a => a.reduce((s, v) => s + v, 0);
  const mean = a => sum(a) / a.length;

  function validateData(c) {
    check(["linear", "gaussian"].includes(c.family), "Poissonian and Gamma were not implemented in the R reference. Choose Linear or Gaussian.");
    check(["constant", "linear", "log-linear"].includes(c.noise), "Choose a supported noise model.");
    check(Number.isInteger(c.seed) && c.seed >= 0 && c.seed <= 2147483647, "Seed must be an integer from 0 to 2147483647.");
    check(Number.isInteger(c.samples) && c.samples >= (c.family === "gaussian" ? 3 : 2) && c.samples <= 10000,
      `Use ${c.family === "gaussian" ? 3 : 2}–10000 x samples; fewer cannot identify this fit.`);
    check(Number.isInteger(c.repeats) && c.repeats >= 1 && c.repeats <= 10000, "Use 1–10000 noise realisations.");
    check(Number.isFinite(c.xmin) && Number.isFinite(c.xmax) && c.xmin >= -10 && c.xmax <= 10 && c.xmin < c.xmax,
      "The sampling range must increase and lie within −10 to 10.");
    check(Number.isFinite(c.noiseMin) && Number.isFinite(c.noiseMax) && c.noiseMin >= -3 && c.noiseMax <= 1 && c.noiseMin <= c.noiseMax,
      "Log10 noise endpoints must be ordered, within −3 to 1.");
    check(c.samples * c.repeats <= 1000000, "Too much stored data: reduce x samples or noise realisations (maximum 1,000,000 observations).");
    return c;
  }

  function createData(config) {
    const c = { ...validateData(config) }, rng = createRng(c.seed);
    const x = grid(c.xmin, c.xmax, c.samples);
    const truth = x.map(t => c.family === "linear" ? t : Astro.normalPdf(t));
    const lo = 10 ** c.noiseMin, hi = 10 ** c.noiseMax;
    const sigma = x.map((_, i) => {
      const t = i / (x.length - 1);
      return c.noise === "constant" ? (lo + hi) / 2
        : c.noise === "linear" ? lo + t * (hi - lo) : 10 ** (c.noiseMin + t * (c.noiseMax - c.noiseMin));
    });
    const observations = Array.from({ length: c.repeats },
      () => truth.map((y, i) => y + sigma[i] * rng.normal()));
    return { config: c, x, truth, sigma, observations };
  }

  function weights(sigma, kind) {
    check(["legacy", "equal", "known"].includes(kind), "Unknown fitting weights.");
    check(sigma.length > 0 && sigma.every(s => Number.isFinite(s) && s > 0), "All observation SDs must be positive and finite.");
    return sigma.map(s => kind === "known" ? 1 / (s * s) : kind === "legacy" ? 1 / s : 1);
  }

  function validateFit(x, y, sigma, w, width, family) {
    check(["linear", "gaussian"].includes(family), "Unsupported fitted model.");
    check(x.length === y.length && x.length === sigma.length && x.length === w.length && x.length >= (family === "gaussian" ? 3 : 2),
      "Fit vectors have incompatible lengths or too few observations.");
    check([x, y, sigma, w].every(a => a.every(Number.isFinite)) && sigma.every(s => s > 0) && w.every(v => v > 0),
      "Fit observations, positive SDs and positive weights must be finite.");
    check(Number.isFinite(width) && width > 0, "The model width / slope must be positive.");
    check(x[0] < x[x.length - 1] && x.every((v, i) => i === 0 || v > x[i - 1]), "x positions must increase strictly.");
  }

  function linearFit(x, y, sigma, w, slope) {
    const sw = sum(w);
    const shift = sum(y.map((v, i) => w[i] * (v - slope * x[i]))) / sw;
    const se = Math.sqrt(sum(w.map((v, i) => (v * sigma[i]) ** 2))) / sw;
    return { shift, amplitude: null, se, amplitudeSE: null, boundary: false, uncertainty: "Exact conditional on slope and fixed weights" };
  }

  // Profile the unrestricted amplitude analytically, then search all resolved
  // centroid peaks. Scaling the template avoids tail underflow.
  function gaussianFit(x, y, sigma, w, width) {
    const lo = x[0], hi = x[x.length - 1];
    function profile(mu) {
      const logs = x.map(v => -0.5 * ((v - mu) / width) ** 2);
      const peak = Math.max(...logs);
      const g = logs.map(v => Math.exp(v - peak));
      let gg = 0, yg = 0;
      for (let i = 0; i < x.length; i++) { gg += w[i] * g[i] * g[i]; yg += w[i] * y[i] * g[i]; }
      return { mu, score: yg * yg / gg, scaledA: yg / gg, peak, g };
    }
    const count = Math.max(49, Math.ceil((hi - lo) / (width * 0.6)) + 1);
    const scan = grid(lo, hi, count).map(profile);
    let best = scan.reduce((a, b) => a.score >= b.score ? a : b);
    const ratio = (Math.sqrt(5) - 1) / 2;
    for (let i = 1; i < count - 1; i++) {
      if (scan[i].score < scan[i - 1].score || scan[i].score < scan[i + 1].score) continue;
      let a = scan[i - 1].mu, b = scan[i + 1].mu;
      let left = profile(b - ratio * (b - a)), right = profile(a + ratio * (b - a));
      for (let k = 0; k < 48; k++) {
        if (left.score > right.score) {
          b = right.mu; right = left; left = profile(b - ratio * (b - a));
        } else {
          a = left.mu; left = right; right = profile(a + ratio * (b - a));
        }
      }
      for (const candidate of [left, right]) if (candidate.score > best.score) best = candidate;
    }
    const amplitude = best.scaledA * Math.exp(-best.peak) * width * Math.sqrt(2 * Math.PI);
    check(Number.isFinite(amplitude), "The Gaussian amplitude is numerically unidentifiable on this x grid. Increase sampling or model width.");
    const fitted = best.g.map(g => best.scaledA * g);
    const boundary = Math.min(best.mu - lo, hi - best.mu) <= 1e-6 * (hi - lo);
    // Implicit differentiation of the score includes the residual-curvature
    // term. J'WJ alone is not the local Hessian for a misspecified width.
    let h00 = 0, h01 = 0, h11 = 0, b00 = 0, b01 = 0, b11 = 0;
    for (let i = 0; i < x.length; i++) {
      const t = (x[i] - best.mu) / (width * width), g = best.g[i];
      const j0 = fitted[i] * t, j1 = g, r = y[i] - fitted[i];
      h00 += w[i] * (j0 * j0 - r * fitted[i] * (t * t - 1 / (width * width)));
      h01 += w[i] * (j0 * j1 - r * g * t);
      h11 += w[i] * j1 * j1;
      const v = (w[i] * sigma[i]) ** 2;
      b00 += v * j0 * j0; b01 += v * j0 * j1; b11 += v * j1 * j1;
    }
    const det = h00 * h11 - h01 * h01;
    let se = null, amplitudeSE = null;
    const regular = !boundary && h00 > 0 && det > 1e-10 * h00 * h11;
    if (regular) {
      const varMu = (h11 * h11 * b00 - 2 * h11 * h01 * b01 + h01 * h01 * b11) / (det * det);
      const varA = (h01 * h01 * b00 - 2 * h00 * h01 * b01 + h00 * h00 * b11) / (det * det);
      if (varMu >= 0 && varA >= 0 && Number.isFinite(varMu) && Number.isFinite(varA)) {
        se = Math.sqrt(varMu);
        amplitudeSE = Math.sqrt(varA) * Math.exp(-best.peak) * width * Math.sqrt(2 * Math.PI);
      }
    }
    return { shift: best.mu, amplitude, se, amplitudeSE, boundary,
      uncertainty: boundary ? "Centroid at search boundary: no symmetric interval"
        : se === null ? "Singular / nonregular fit: uncertainty unavailable" : "Local noise propagation, conditional on width and fixed weights" };
  }

  function fit(x, y, sigma, width, family, w) {
    validateFit(x, y, sigma, w, width, family);
    return family === "linear" ? linearFit(x, y, sigma, w, width) : gaussianFit(x, y, sigma, w, width);
  }

  function predict(x, fitted, width, family) {
    return x.map(t => family === "linear" ? width * t + fitted.shift : fitted.amplitude * Astro.normalPdf(t, fitted.shift, width));
  }

  function residualWeights(x, y, initial, width, family, base) {
    const f = predict(x, initial, width, family);
    const residual = y.map((v, i) => Math.abs(v - f[i]));
    const fm = mean(f), rm = mean(residual);
    const ss = sum(f.map(v => (v - fm) ** 2));
    check(ss > 1e-20, "Residual reweighting cannot regress on a constant fitted response.");
    const slope = sum(f.map((v, i) => (v - fm) * (residual[i] - rm))) / ss;
    const predicted = f.map(v => rm + slope * (v - fm));
    check(predicted.every(v => Number.isFinite(v) && v > 1e-12),
      "Residual reweighting predicts a zero or negative noise scale. No weights were clipped: turn it off or use known inverse-variance weights.");
    const result = base.map((v, i) => v / (predicted[i] ** 2));
    const scale = Math.max(...result);
    return result.map(v => v / scale);
  }

  function summarize(fits, reference) {
    const values = fits.map(f => f.shift), m = mean(values);
    const sd = values.length > 1 ? Math.sqrt(sum(values.map(v => (v - m) ** 2)) / (values.length - 1)) : null;
    const intervals = fits.filter(f => f.se !== null);
    const coverage = centre => intervals.length
      ? intervals.filter(f => Math.abs(f.shift - centre) <= Z95 * f.se).length / intervals.length : null;
    return { mean: m, sd, mcse: sd === null ? null : sd / Math.sqrt(fits.length),
      reference: reference.shift, referenceSE: reference.se, coverageTruth: coverage(0),
      coverageReference: coverage(reference.shift), intervalCount: intervals.length,
      boundaryCount: fits.filter(f => f.boundary).length };
  }

  function compare(data, settings) {
    const { widthMin, widthMax, models, weighting, residual } = settings;
    check(Number.isInteger(models) && models >= 2 && models <= 8, "Use 2–8 fitted models.");
    check(Number.isFinite(widthMin) && Number.isFinite(widthMax) && widthMin >= 0.1 && widthMax <= 1.9 && widthMin <= widthMax,
      "Model slope / width endpoints must be ordered, within 0.1–1.9.");
    check(typeof residual === "boolean", "Residual reweighting must be true or false.");
    const widths = grid(widthMin, widthMax, models), family = data.config.family;
    const cost = data.x.length * data.observations.length * models *
      (family === "gaussian" ? Math.max(49, Math.ceil((data.config.xmax - data.config.xmin) / (widthMin * 0.6)) + 1) + 100 : 1) * (residual ? 2 : 1);
    check(cost <= 60000000, "This combination exceeds the browser fitting budget. Reduce samples, realisations, models, or widen the Gaussian template (no data were thinned).");
    const base = weights(data.sigma, weighting);
    const run = (y, width, reweight) => {
      const first = fit(data.x, y, data.sigma, width, family, base);
      if (!reweight) return first;
      const w = residualWeights(data.x, y, first, width, family, base);
      const second = fit(data.x, y, data.sigma, width, family, w);
      return { ...second, se: null, amplitudeSE: null,
        uncertainty: "Data-dependent heuristic weights: no calibrated analytic interval" };
    };
    return widths.map(width => {
      // The reference is always the noiseless, base-weight fit: fitting an
      // absolute-residual noise model to noiseless residuals is not meaningful.
      const reference = run(data.truth, width, false);
      const fits = data.observations.map(y => run(y, width, residual));
      return { width, reference, fits, summary: summarize(fits, reference) };
    });
  }

  function histogram(values, bins = 36) {
    check(values.length > 0 && values.every(Number.isFinite), "Histogram needs finite fitted shifts.");
    const [lo, hi] = values.reduce(([a, b], v) => [Math.min(a, v), Math.max(b, v)], [-0.2, 0.2]);
    const step = (hi - lo) / bins;
    const counts = Array(bins).fill(0);
    for (const value of values) counts[Math.min(bins - 1, Math.floor((value - lo) / step))]++;
    return { lo, hi, step, points: counts.map((n, i) => [lo + (i + 0.5) * step, n / (values.length * step)]) };
  }

  // Paired groups share each realisation. Regress group means and propagate
  // whole-realisation slopes, rather than pretending all groups independent.
  function biasTrend(groups) {
    check(groups.length >= 2 && groups[0].fits.length >= 2, "Trend bands need two models and two realisations.");
    const widths = groups.map(g => g.width), xm = mean(widths), ss = sum(widths.map(v => (v - xm) ** 2));
    check(ss > 0, "Trend bands require distinct model parameters.");
    const rows = groups[0].fits.map((_, j) => {
      const values = groups.map(g => g.fits[j].shift);
      const interceptAtMean = mean(values);
      const slope = sum(values.map((v, i) => (widths[i] - xm) * v)) / ss;
      return { interceptAtMean, slope };
    });
    const x = grid(Math.min(...widths), Math.max(...widths), 120);
    return x.map(v => {
      const predictions = rows.map(r => r.interceptAtMean + r.slope * (v - xm));
      const m = mean(predictions), sd = Math.sqrt(sum(predictions.map(p => (p - m) ** 2)) / (rows.length - 1));
      return { x: v, mean: m, meanSE: sd / Math.sqrt(rows.length), sd };
    });
  }

  return { Z95, grid, validateData, createData, weights, fit, predict, residualWeights, summarize, compare, histogram, biasTrend };
});
