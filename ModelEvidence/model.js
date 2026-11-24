(function (root) {
  "use strict";
  const A = typeof module === "object" && module.exports ? require("../TeachingTools/statistics.js") : root.Astro;
  const defaults = { amplitude: 0.5, noise: 1, prior: 2, seed: 42 };
  function lineEvidence(estimate, error, priorSd) {
    A.check(Number.isFinite(estimate) && Number.isFinite(error) && error > 0 &&
      Number.isFinite(priorSd) && priorSd > 0, "Estimate must be finite; measurement and proper prior SDs must be positive.");
    const ratio = (priorSd / error) ** 2, z2 = (estimate / error) ** 2;
    const logNull = -Math.log(error) - 0.5 * Math.log(2 * Math.PI) - z2 / 2;
    const logAlternative = -0.5 * Math.log(2 * Math.PI * (error ** 2 + priorSd ** 2)) -
      estimate ** 2 / (2 * (error ** 2 + priorSd ** 2));
    const logBf = -0.5 * Math.log1p(ratio) + 0.5 * z2 * (ratio / (1 + ratio));
    A.check([logNull, logAlternative, logBf].every(Number.isFinite), "The evidence exceeds numerical precision.");
    const probability = logBf >= 0 ? 1 / (1 + Math.exp(-logBf)) : Math.exp(logBf) / (1 + Math.exp(logBf));
    return { logBf, probability, logMaximumRatio: z2 / 2, logNull, logAlternative,
      posteriorMean: estimate * ratio / (1 + ratio), posteriorSd: error * Math.sqrt(ratio / (1 + ratio)) };
  }
  function compute(parameters = {}) {
    const p = { ...defaults, ...parameters };
    A.check(Number.isFinite(p.amplitude) && p.amplitude >= -2 && p.amplitude <= 3, "True amplitude must be between -2 and 3.");
    A.check(Number.isFinite(p.noise) && p.noise >= 0.3 && p.noise <= 2, "Flux uncertainty must be between 0.3 and 2.");
    A.check(Number.isFinite(p.prior) && p.prior >= 0.1 && p.prior <= 10, "The proper amplitude-prior SD must be between 0.1 and 10.");
    const rng = A.createRng(p.seed), n = 70;
    const wavelength = Array.from({ length: n }, (_, i) => 500 + 10 * i / (n - 1));
    const template = wavelength.map(x => Math.exp(-0.5 * ((x - 505) / 0.35) ** 2));
    const flux = template.map(g => 1 + p.amplitude * g + p.noise * rng.normal());
    const sumSquares = template.reduce((sum, g) => sum + g * g, 0);
    const estimate = template.reduce((sum, g, i) => sum + g * (flux[i] - 1), 0) / sumSquares;
    const error = p.noise / Math.sqrt(sumSquares);
    const evidence = lineEvidence(estimate, error, p.prior);
    const residualSquares = flux.reduce((sum, value, i) => sum + (value - 1 - estimate * template[i]) ** 2, 0);
    const commonLogFactor = -n / 2 * Math.log(2 * Math.PI * p.noise ** 2) -
      residualSquares / (2 * p.noise ** 2) + 0.5 * Math.log(2 * Math.PI * error ** 2);
    return { wavelength, template, flux, estimate, error, evidence, residualSquares,
      logEvidenceNull: commonLogFactor + evidence.logNull,
      logEvidenceLine: commonLogFactor + evidence.logAlternative };
  }
  const api = { defaults, lineEvidence, compute };
  if (typeof module === "object" && module.exports) module.exports = api;
  else root.ModelEvidenceModel = api;
})(globalThis);
