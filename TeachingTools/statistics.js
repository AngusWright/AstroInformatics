(function (root) {
  "use strict";
  class ValidationError extends RangeError {}
  function check(condition, message) { if (!condition) throw new ValidationError(message); }
  const LOG_2PI = Math.log(2 * Math.PI);
  const EPS = 4e-15;

  function createRng(seed) {
    check(Number.isInteger(seed) && seed >= 0 && seed <= 2147483647,
      "Seed must be an integer between 0 and 2147483647.");
    let state = seed >>> 0, spare = null;
    function uniform() {
      state = (state + 0x6d2b79f5) >>> 0;
      let t = state;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return (((t ^ (t >>> 14)) >>> 0) + 0.5) / 4294967296;
    }
    function normal() {
      if (spare !== null) { const result = spare; spare = null; return result; }
      const radius = Math.sqrt(-2 * Math.log(uniform())), angle = 2 * Math.PI * uniform();
      spare = radius * Math.sin(angle);
      return radius * Math.cos(angle);
    }
    function gamma(shape, rate) {
      check(Number.isFinite(shape) && shape > 0 && Number.isFinite(rate) && rate > 0,
        "Gamma shape and rate must be positive and finite.");
      if (shape < 1) return gamma(shape + 1, rate) * uniform() ** (1 / shape);
      const d = shape - 1 / 3, c = 1 / Math.sqrt(9 * d);
      for (;;) {
        const x = normal(), base = 1 + c * x;
        if (base <= 0) continue;
        const v = base ** 3, u = uniform();
        if (u < 1 - 0.0331 * x ** 4 ||
            Math.log(u) < x * x / 2 + d * (1 - v + Math.log(v))) return d * v / rate;
      }
    }
    return { uniform, normal, gamma };
  }

  function logGamma(z) {
    check(Number.isFinite(z) && z > 0, "Gamma argument must be positive and finite.");
    const coefficients = [
      676.5203681218851, -1259.1392167224028, 771.3234287776531,
      -176.6150291621406, 12.507343278686905, -0.13857109526572012,
      9.984369578019572e-6, 1.5056327351493116e-7
    ];
    if (z < 0.5) return Math.log(Math.PI) - Math.log(Math.sin(Math.PI * z)) - logGamma(1 - z);
    z -= 1;
    let sum = 0.9999999999998099;
    coefficients.forEach((value, i) => { sum += value / (z + i + 1); });
    const t = z + 7.5;
    return LOG_2PI / 2 + (z + 0.5) * Math.log(t) - t + Math.log(sum);
  }

  function gammaTails(a, x, logNormaliser = logGamma(a)) {
    if (x <= 0) return [0, 1];
    if (x === Infinity) return [1, 0];
    const factor = Math.exp(a * Math.log(x) - x - logNormaliser);
    if (x < a + 1) {
      let term = 1 / a, sum = term;
      for (let i = 1; i <= 10000; i++) {
        term *= x / (a + i); sum += term;
        if (term <= sum * EPS) {
          const p = Math.min(1, factor * sum);
          return [p, 1 - p];
        }
      }
    } else {
      let b = x + 1 - a, c = 1e300, d = 1 / b, h = d;
      for (let i = 1; i <= 10000; i++) {
        const an = -i * (i - a);
        b += 2; d = an * d + b; c = b + an / c;
        if (Math.abs(d) < 1e-300) d = d < 0 ? -1e-300 : 1e-300;
        if (Math.abs(c) < 1e-300) c = c < 0 ? -1e-300 : 1e-300;
        d = 1 / d;
        const delta = d * c;
        h *= delta;
        if (Math.abs(delta - 1) < EPS) {
          const q = Math.min(1, factor * h);
          return [1 - q, q];
        }
      }
    }
    throw new ValidationError("The Gamma CDF did not converge for these parameters.");
  }

  function betaFraction(a, b, x) {
    const floor = value => Math.abs(value) < 1e-300 ? (value < 0 ? -1e-300 : 1e-300) : value;
    let c = 1, d = 1 / floor(1 - (a + b) * x / (a + 1)), h = d;
    for (let m = 1; m <= 10000; m++) {
      const m2 = 2 * m;
      let aa = m * (b - m) * x / ((a + m2 - 1) * (a + m2));
      d = 1 / floor(1 + aa * d); c = floor(1 + aa / c); h *= d * c;
      aa = -(a + m) * (a + b + m) * x / ((a + m2) * (a + m2 + 1));
      d = 1 / floor(1 + aa * d); c = floor(1 + aa / c);
      const delta = d * c; h *= delta;
      if (Math.abs(delta - 1) < EPS) return h;
    }
    throw new ValidationError("The Beta CDF did not converge for these parameters.");
  }

  function betaTails(a, b, x, logBeta) {
    if (x <= 0) return [0, 1];
    if (x >= 1) return [1, 0];
    const factor = Math.exp(a * Math.log(x) + b * Math.log1p(-x) - logBeta);
    if (x < (a + 1) / (a + b + 2)) {
      const p = Math.min(1, factor * betaFraction(a, b, x) / a);
      return [p, 1 - p];
    }
    const q = Math.min(1, factor * betaFraction(b, a, 1 - x) / b);
    return [1 - q, q];
  }

  function normalPdf(x, mean = 0, sd = 1) {
    return Math.exp(-0.5 * ((x - mean) / sd) ** 2 - Math.log(sd) - LOG_2PI / 2);
  }

  function normalCdf(x, mean = 0, sd = 1) {
    const z = (x - mean) / sd;
    const tail = gammaTails(0.5, z * z / 2)[1] / 2;
    return z < 0 ? tail : 1 - tail;
  }

  function normalQuantile(p, mean = 0, sd = 1) {
    check(Number.isFinite(p) && p >= 0 && p <= 1, "Probability must be between zero and one.");
    if (p === 0) return -Infinity;
    if (p === 1) return Infinity;
    const a = [-39.69683028665376, 220.9460984245205, -275.9285104469687, 138.357751867269,
      -30.66479806614716, 2.506628277459239];
    const b = [-54.47609879822406, 161.5858368580409, -155.6989798598866, 66.80131188771972, -13.28068155288572];
    const c = [-0.007784894002430293, -0.3223964580411365, -2.400758277161838,
      -2.549732539343734, 4.374664141464968, 2.938163982698783];
    const d = [0.007784695709041462, 0.3224671290700398, 2.445134137142996, 3.754408661907416];
    const poly = (coefficients, x) => coefficients.reduce((sum, coefficient) => sum * x + coefficient, 0);
    let z;
    if (p < 0.02425 || p > 0.97575) {
      const q = Math.sqrt(-2 * (p < 0.5 ? Math.log(p) : Math.log1p(-p)));
      z = poly(c, q) / (poly(d, q) * q + 1);
      if (p > 0.5) z = -z;
    } else {
      const q = p - 0.5, r = q * q;
      z = q * poly(a, r) / (poly(b, r) * r + 1);
    }
    z -= (normalCdf(z) - p) / normalPdf(z);
    return mean + sd * z;
  }

  function invert(p, tails, lo, hi) {
    check(Number.isFinite(p) && p > 0 && p < 1, "Interior quantiles need 0 < probability < 1.");
    const above = x => p < 0.5 ? tails(x)[0] >= p : tails(x)[1] <= 1 - p;
    while (!above(hi)) {
      hi *= 2;
      check(Number.isFinite(hi), "The quantile exceeds numerical precision.");
    }
    for (let i = 0; i < 160; i++) {
      const mid = lo + (hi - lo) / 2;
      if (mid === lo || mid === hi) return mid;
      if (above(mid)) hi = mid; else lo = mid;
    }
    const result = lo + (hi - lo) / 2;
    check(Math.abs(tails(result)[p < 0.5 ? 0 : 1] - Math.min(p, 1 - p)) <= 1e-8 * Math.min(p, 1 - p),
      "The quantile cannot be resolved at floating-point precision.");
    return result;
  }

  function normal(mean, sd) {
    check(Number.isFinite(mean) && Number.isFinite(sd) && sd > 0, "Normal mean must be finite and SD positive.");
    return {
      label: `Normal: mean ${mean.toPrecision(4)}, SD ${sd.toPrecision(4)}`, mean, sd,
      range: [mean - 4.5 * sd, mean + 4.5 * sd], focus: [{ location: mean, scale: sd }],
      density: x => normalPdf(x, mean, sd), cdf: x => normalCdf(x, mean, sd),
      quantile: p => normalQuantile(p, mean, sd), draw: rng => mean + sd * rng.normal()
    };
  }

  function beta(a, b) {
    check(Number.isFinite(a) && a > 0 && Number.isFinite(b) && b > 0, "Beta shapes must be positive and finite.");
    const logBeta = logGamma(a) + logGamma(b) - logGamma(a + b);
    const tails = x => betaTails(a, b, x, logBeta);
    return {
      label: `Beta(${a.toPrecision(4)}, ${b.toPrecision(4)})`,
      mean: a / (a + b), sd: Math.sqrt(a * b / ((a + b) ** 2 * (a + b + 1))),
      density(x) {
        if (x < 0 || x > 1) return 0;
        if (x === 0) return a < 1 ? Infinity : a === 1 ? b : 0;
        if (x === 1) return b < 1 ? Infinity : b === 1 ? a : 0;
        return Math.exp((a - 1) * Math.log(x) + (b - 1) * Math.log1p(-x) - logBeta);
      },
      cdf: x => tails(x)[0],
      quantile(p) {
        if (p === 0 || p === 1) return p;
        return invert(p, tails, 0, 1);
      },
      draw(rng) { const x = rng.gamma(a, 1), y = rng.gamma(b, 1); return x / (x + y); }
    };
  }

  function gamma(shape, rate) {
    check(Number.isFinite(shape) && shape > 0 && Number.isFinite(rate) && rate > 0,
      "Gamma shape and rate must be positive and finite.");
    const mean = shape / rate, sd = Math.sqrt(shape) / rate;
    const logNormaliser = logGamma(shape), constant = shape * Math.log(rate) - logNormaliser;
    const tails = x => gammaTails(shape, x, logNormaliser);
    return {
      label: `Gamma: shape ${shape.toPrecision(4)}, rate ${rate.toPrecision(4)}`,
      mean, sd, range: [0, mean + 7 * sd], focus: [{ location: Math.max(0, (shape - 1) / rate), scale: sd }],
      density(x) {
        if (x < 0 || x === Infinity) return 0;
        if (x === 0) return shape === 1 ? rate : shape < 1 ? Infinity : 0;
        return Math.exp(constant + (shape - 1) * Math.log(x) - rate * x);
      },
      cdf: x => tails(x * rate)[0],
      quantile(p) {
        if (p === 0) return 0;
        if (p === 1) return Infinity;
        return invert(p, tails, 0, Math.max(1, shape + 10 * Math.sqrt(shape))) / rate;
      },
      draw: rng => rng.gamma(shape, rate)
    };
  }

  function summarize(values) {
    check(values.length >= 2 && values.every(Number.isFinite), "At least two finite values are required.");
    let mean = 0, sum = 0;
    values.forEach((x, i) => { const delta = x - mean; mean += delta / (i + 1); sum += delta * (x - mean); });
    return { mean, sd: Math.sqrt(sum / (values.length - 1)) };
  }

  function acf(values, maxLag = 40) {
    const { mean, sd } = summarize(values);
    check(sd > 0, "Autocorrelation is undefined without variation.");
    check(Number.isInteger(maxLag) && maxLag >= 0, "Maximum lag must be a non-negative integer.");
    const denominator = values.reduce((sum, x) => sum + (x - mean) ** 2, 0);
    return Array.from({ length: Math.min(maxLag, values.length - 1) + 1 }, (_, lag) =>
      values.slice(lag).reduce((sum, x, i) => sum + (x - mean) * (values[i] - mean), 0) / denominator);
  }

  function quantile(values, probability) {
    check(values.length > 0 && values.every(Number.isFinite), "Quantiles need finite observations.");
    check(Number.isFinite(probability) && probability >= 0 && probability <= 1, "Probability must be between zero and one.");
    const sorted = values.slice().sort((a, b) => a - b), position = (sorted.length - 1) * probability;
    const lower = Math.floor(position), upper = Math.ceil(position);
    return sorted[lower] + (position - lower) * (sorted[upper] - sorted[lower]);
  }

  function histogram(values, breaks) {
    check(values.every(Number.isFinite), "Histogram observations must be finite.");
    check(breaks.length >= 2 && breaks.every(Number.isFinite) &&
      breaks.every((value, i) => i === 0 || value > breaks[i - 1]), "Histogram edges must be finite and increasing.");
    const counts = Array(breaks.length - 1).fill(0);
    let excluded = 0;
    for (const value of values) {
      if (value < breaks[0] || value > breaks[breaks.length - 1]) { excluded++; continue; }
      let lo = 0, hi = breaks.length - 1;
      while (hi - lo > 1) {
        const mid = Math.floor((lo + hi) / 2);
        if (value < breaks[mid]) hi = mid;
        else lo = mid;
      }
      counts[Math.min(lo, counts.length - 1)]++;
    }
    return {
      counts, excluded,
      mids: counts.map((_, i) => (breaks[i] + breaks[i + 1]) / 2),
      density: counts.map((count, i) => values.length ? count / (values.length * (breaks[i + 1] - breaks[i])) : 0),
      widths: counts.map((_, i) => breaks[i + 1] - breaks[i])
    };
  }

  const api = { ValidationError, check, createRng, logGamma, normalPdf, normalCdf, normalQuantile,
    normal, beta, gamma, summarize, acf, quantile, histogram };
  if (typeof module !== "undefined" && module.exports) module.exports = api;
  root.Astro = api;
})(globalThis);
