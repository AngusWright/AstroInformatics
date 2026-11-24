(function () {
  "use strict";
  const U = AstroUI, M = DensityEstimationModel;
  let cached = null;
  const presets = { under: { bandwidth: 0.05, width: 0.1 }, over: { bandwidth: 1, width: 1.5 },
    middle: { bandwidth: 0.25, width: 0.4 } };
  window.lesson = AstroLesson.mount({
    presets,
    compute(p) {
      M.validate(p);
      if (!cached || cached.seed !== p.seed) cached = M.prepare(p.seed);
      return M.compute(p, cached);
    },
    render(result, p) {
      U.$("objects").textContent = result.values.length;
      U.$("peaks").textContent = result.peaks.length;
      U.$("hist-area").textContent = result.histogramArea.toFixed(3);
      U.$("outside-mass").textContent = (100 * result.outsideMass).toPrecision(3);
      const truth = result.x.map(M.populationDensity);
      const yMax = Math.max(...result.histogram.density, ...result.density, ...truth) * 1.12;
      const plot = U.plot(U.$("density-plot"), result.limits, [0, yMax], "Stellar colour (illustrative units)", "Probability density");
      plot.clip();
      U.bars(plot, result.histogram.mids, result.histogram.density, result.histogram.widths, "#cbd8eb");
      U.line(plot, result.x.map((x, i) => [x, truth[i]]), "#087f8c", 3, [7, 5]);
      U.line(plot, result.x.map((x, i) => [x, result.density[i]]), "#5059a6", 3);
      result.values.forEach(value => U.line(plot, [[value, 0], [value, yMax * 0.045]], "#d5653477", 1));
      plot.ctx.restore();
      U.$("data-note").textContent = `${p.n} identical observations for every smoothing choice. Bin width ${p.width}, origin ${p.origin} × width; Gaussian bandwidth ${p.bandwidth}. All histogram bins are in view; KDE tails are not renormalised.`;
    }
  });
})();
