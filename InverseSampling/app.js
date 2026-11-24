(function () {
  "use strict";
  const A = Astro, U = AstroUI, M = InverseSamplingModel;
  let cached = null;
  const presets = { loguniform: { family: "power", slope: 1 }, steep: { family: "power", slope: 3 } };
  window.lesson = AstroLesson.mount({
    steps: 100, presets,
    compute(p) {
      M.validate(p);
      if (!cached || cached.seed !== p.seed) cached = M.prepare(p.seed);
      return M.compute(p, cached);
    },
    render(result, p, frame) {
      const current = M.atFrame(result, frame), target = result.target;
      const z = U.linspace(...result.limits, 700), cdf = z.map(value => target.cdf(10 ** value));
      const mass = Math.max(0, Math.min(1, cdf[cdf.length - 1] - cdf[0]));
      U.$("revealed").textContent = current.count;
      U.$("ks").textContent = current.discrepancy === null ? "No draws yet" : U.number(current.discrepancy);
      U.$("target-mass").textContent = mass.toFixed(6);
      const upper = U.plot(U.$("cdf-plot"), result.limits, [0, 1.03], "log10 target value", "Cumulative probability");
      upper.clip();
      U.line(upper, z.map((value, i) => [value, cdf[i]]), "#5059a6", 3);
      for (let i = Math.max(0, current.count - 12); i < current.count; i++) {
        const x = Math.log10(current.values[i]), u = current.uniforms[i];
        U.line(upper, [[result.limits[0], u], [x, u], [x, 0]], "#d5653488", 1.5);
        U.points(upper, [[x, u]], "#087f8c", 3.2);
      }
      upper.ctx.restore();
      const histogram = A.histogram(current.values.map(Math.log10), U.linspace(...result.limits, 45));
      const density = z.map(value => M.logDensity(target, value));
      const lower = U.plot(U.$("density-plot"), result.limits, [0, Math.max(...histogram.density, ...density) * 1.12],
        "log10 sampled value", "Density per log10 unit");
      lower.clip();
      U.bars(lower, histogram.mids, histogram.density, histogram.widths, "#e9ab8e");
      U.line(lower, z.map((value, i) => [value, density[i]]), "#087f8c", 2.5);
      lower.ctx.restore();
      U.$("target-description").textContent = `${target.label}. Read a guide horizontally from U to the curve, then down to log10(x).`;
      U.$("tail-note").textContent = `Histogram denominator: all ${current.count} revealed draws; ${histogram.excluded} outside view. Target probability outside view: ${(100 * (1 - mass)).toPrecision(3)}%. No renormalisation.`;
    }
  });
})();
