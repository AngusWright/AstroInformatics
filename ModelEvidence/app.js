(function () {
  "use strict";
  const M = ModelEvidenceModel, U = AstroUI, A = Astro, ln10 = Math.log(10);
  window.lesson = AstroLesson.mount({
    compute: M.compute,
    presets: { "noise-only": { amplitude: 0 }, "clear-line": { amplitude: 2 }, broad: { prior: 10 } },
    render(result, p) {
      const e = result.evidence;
      U.$("estimate").textContent = `${result.estimate.toFixed(2)} ± ${result.error.toFixed(2)}`;
      U.$("log-bf").textContent = (e.logBf / ln10).toFixed(2);
      U.$("probability").textContent = `${(100 * e.probability).toFixed(1)}%`;
      U.$("log-ml").textContent = (e.logMaximumRatio / ln10).toFixed(2);
      U.$("comparison-note").textContent = `H₁ uses the proper prior Normal(0, ${p.prior.toFixed(1)}²). ` +
        `The log₁₀ predictive densities at the measured amplitude are ${(e.logNull / ln10).toFixed(3)} (H₀) and ` +
        `${(e.logAlternative / ln10).toFixed(3)} (H₁). Their difference is the log₁₀ Bayes factor.`;
      U.table("evidence-table", [
        ["log₁₀ p(full spectrum | H₀)", (result.logEvidenceNull / ln10).toFixed(6)],
        ["log₁₀ p(full spectrum | H₁)", (result.logEvidenceLine / ln10).toFixed(6)],
        ["log₁₀ p(Â | H₀)", (e.logNull / ln10).toFixed(6)],
        ["log₁₀ p(Â | H₁)", (e.logAlternative / ln10).toFixed(6)],
        ["A | data, H₁: mean ± SD", `${e.posteriorMean.toFixed(4)} ± ${e.posteriorSd.toFixed(4)}`],
        ["Prior model odds H₁ : H₀", "1 : 1"]
      ]);
      const fitted = result.template.map(g => 1 + result.estimate * g);
      const spectrum = U.plot(U.$("spectrum-plot"), [500, 510], U.extent([...result.flux, ...fitted, 1]),
        "Wavelength (nm)", "Flux");
      spectrum.clip();
      U.line(spectrum, [[500, 1], [510, 1]], "#087f8c", 2, [6, 4]);
      U.line(spectrum, result.wavelength.map((x, i) => [x, fitted[i]]), "#5059a6", 2.5);
      U.points(spectrum, result.wavelength.map((x, i) => [x, result.flux[i]]), "#d56534", 3);
      spectrum.ctx.restore();
      const alternativeSd = Math.hypot(result.error, p.prior);
      const half = Math.max(3.5 * alternativeSd, Math.abs(result.estimate) * 1.15);
      const xs = [...U.linspace(-half, half, 700), ...U.linspace(-5 * result.error, 5 * result.error, 400), result.estimate]
        .sort((a, b) => a - b).filter(x => x >= -half && x <= half);
      const peak = A.normalPdf(0, 0, result.error);
      const predictive = U.plot(U.$("predictive-plot"), [-half, half], [0, peak * 1.08],
        "Possible measured amplitude Â", "Predictive density");
      predictive.clip();
      U.line(predictive, xs.map(x => [x, A.normalPdf(x, 0, result.error)]), "#087f8c", 2.5);
      U.line(predictive, xs.map(x => [x, A.normalPdf(x, 0, alternativeSd)]), "#5059a6", 2.5);
      U.line(predictive, [[result.estimate, 0], [result.estimate, peak * 1.08]], "#d56534", 2, [6, 4]);
      U.points(predictive, [[result.estimate, Math.exp(e.logNull)]], "#087f8c", 4);
      U.points(predictive, [[result.estimate, Math.exp(e.logAlternative)]], "#5059a6", 4);
      predictive.ctx.restore();
    }
  });
})();
