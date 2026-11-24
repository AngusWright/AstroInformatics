(function () {
  "use strict";
  const U = AstroUI, M = CorrelationModel;
  let base = null;
  const formatR = value => value === null ? "Undefined" : value.toFixed(2);
  function colour(z, alpha) {
    const stops = [[8, 127, 140], [182, 221, 210], [235, 154, 84]];
    const half = Math.min(1, Math.floor(z * 2)), t = z * 2 - half;
    const rgb = stops[half].map((value, i) => Math.round(value + t * (stops[half + 1][i] - value)));
    return `rgba(${rgb.join(",")},${alpha})`;
  }
  function population(frame, result, x, y) {
    const bins = Array.from({ length: 60 }, () => ({ all: [], selected: [] }));
    result.distance.forEach((z, i) => {
      const bucket = bins[Math.min(59, Math.floor(60 * z))], point = [x[i], y[i]];
      bucket.all.push(point);
      if (result.selected[i]) bucket.selected.push(point);
    });
    bins.forEach((bucket, i) => U.points(frame, bucket.all, colour((i + 0.5) / 60, 0.18), 1.8));
    bins.forEach((bucket, i) => U.points(frame, bucket.selected, colour((i + 0.5) / 60, 0.9), 2.8));
  }
  function regression(frame, stats, range, colour, dash = []) {
    if (stats.slope === null) return;
    U.line(frame, range.map(x => [x, stats.intercept + stats.slope * x]), colour, 2, dash);
  }
  function render(result) {
    U.$("full-r").textContent = formatR(result.full.correlation);
    U.$("slice-r").textContent = result.slice.n < 3 ? "Too few objects" : formatR(result.slice.correlation);
    U.$("intrinsic-r").textContent = formatR(result.intrinsic.correlation);
    U.$("selected-count").textContent = `${result.slice.n} / ${M.SIZE}`;
    const xRange = U.extent(result.x), yRange = U.extent(result.y);
    const flux = U.plot(U.$("flux-plot"), xRange, yRange, "Log flux in band X", "Log flux in band Y");
    flux.clip(); population(flux, result, result.x, result.y);
    regression(flux, result.full, xRange, "#7f8a9b", [6, 4]);
    regression(flux, result.slice, xRange, "#5059a6");
    flux.ctx.restore();
    const residualRange = U.extent([...result.intrinsicX, ...result.intrinsicY]);
    const intrinsic = U.plot(U.$("intrinsic-plot"), residualRange, residualRange, "Intrinsic X = X + 2Z", "Intrinsic Y = Y + 2Z");
    intrinsic.clip(); population(intrinsic, result, result.intrinsicX, result.intrinsicY); intrinsic.ctx.restore();
    U.$("population-note").textContent = `Known-model population correlations: whole sample ${formatR(result.theoryFull)}; ` +
      `finite distance slice ${formatR(result.theorySlice)} (Z from ${result.lower.toFixed(3)} to ${result.upper.toFixed(3)}). ` +
      "The cards show sample correlations, not these population expectations.";
    U.$("variance-note").textContent = result.intrinsic.correlation === null ?
      "No intrinsic variation: Pearson r is undefined when either variable has zero variance. Exact-distance conditional variables are constants here." :
      "Independent intrinsic terms have population correlation zero. Nonzero sample r and residual distance variation in a finite slice are expected.";
  }
  globalThis.lesson = AstroLesson.mount({
    compute(parameters) {
      M.validate(parameters);
      if (base === null || base.seed !== parameters.seed) base = M.prepare(parameters.seed);
      return M.compute(parameters, base);
    },
    render,
    presets: {
      whole: { centre: 0.5, width: 1 },
      narrow: { centre: 0.5, width: 0.05 },
      deterministic: { scatter: 0, centre: 0.5, width: 0.2 }
    }
  });
})();
