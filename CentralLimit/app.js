(function () {
  "use strict";
  const A = Astro, U = AstroUI, M = CentralLimitModel;
  let cached = null;
  const presets = {
    one: { n: 1 },
    large: { n: 100 },
    cauchy: { population: "cauchy", n: 100 }
  };
  window.lesson = AstroLesson.mount({
    steps: 100, presets,
    compute(p) {
      M.validate(p);
      if (!cached || cached.seed !== p.seed || cached.population !== p.population) cached = M.prepare(p.seed, p.population);
      return M.compute(p, cached);
    },
    render(result, p, frame) {
      const current = M.atFrame(result, frame), population = result.population;
      U.$("stack-n").textContent = current.n;
      U.$("experiments").textContent = p.repetitions;
      U.$("mean-sd").textContent = current.finite ? U.number(current.meanSd) : "Undefined";
      U.$("empirical-iqr").textContent = U.number(current.empiricalIqr);
      U.$("moment-note").textContent = current.finite ?
        `Population mean ${U.number(population.mean)}, variance ${U.number(population.variance)}. Raw averages narrow as 1/√n; the lower figure standardises that narrowing away to test the CLT shape.` :
        "Cauchy population mean and variance are undefined. Every raw average is still exactly Cauchy: its central 50% width is 2 at every stack size. A finite-run IQR is not a population moment.";
      const populationX = U.linspace(...population.limits, 600);
      const populationY = populationX.map(x => M.populationDensity(p.population, x));
      const populationHistogram = A.histogram(result.data.observations, U.linspace(...population.limits, 65));
      const upper = U.plot(U.$("population-plot"), population.limits,
        [0, Math.max(...populationHistogram.density, ...populationY) * 1.12], "Individual observation", "Density");
      upper.clip();
      U.bars(upper, populationHistogram.mids, populationHistogram.density, populationHistogram.widths, "#cbd8eb");
      U.line(upper, populationX.map((x, i) => [x, populationY[i]]), "#087f8c", 2.5);
      upper.ctx.restore();
      const limits = current.finite ? [-4, 6] : [-10, 10], x = U.linspace(...limits, 600);
      const reference = x.map(value => current.finite ? A.normalPdf(value) : M.populationDensity("cauchy", value));
      const histogram = A.histogram(current.statistic, U.linspace(...limits, 60));
      const lower = U.plot(U.$("sampling-plot"), limits, [0, Math.max(...histogram.density, ...reference) * 1.12],
        current.finite ? "Standardised average Z" : "Raw average (Cauchy)", "Density");
      lower.clip();
      U.bars(lower, histogram.mids, histogram.density, histogram.widths, "#e9ab8e");
      U.line(lower, x.map((value, i) => [value, reference[i]]), "#5059a6", 2.5);
      lower.ctx.restore();
      U.$("sampling-title").textContent = `Independent averages: n = ${current.n}`;
      U.$("reference-label").textContent = current.finite ? "Purple: Normal(0,1) CLT reference" : "Purple: exact standard Cauchy reference";
      U.$("population-tails").textContent = `${population.label}. Outside this view: ${populationHistogram.excluded} / 8000 observations; no tail renormalisation.`;
      U.$("sampling-tails").textContent = current.finite ?
        `Z = √n (average − μ) / σ. Outside this view: ${histogram.excluded} / ${p.repetitions} experiments. The Gaussian curve is an asymptotic reference.` :
        `No Gaussian standardisation is possible. Outside this view: ${histogram.excluded} / ${p.repetitions} raw averages. The reference does not shrink with n.`;
    }
  });
})();
