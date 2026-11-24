(function () {
  "use strict";
  const U = AstroUI, M = PosteriorPredictiveModel;
  const names = ["mean", "star", "measurement"];
  const labels = ["Population mean", "Next true stellar velocity", "Next measured velocity"];
  const colours = ["#5059A6", "#087F8C", "#D56534"];
  function render(result, p) {
    const posterior = result.posterior, centre = posterior.mean;
    const sds = names.map(name => Math.sqrt(posterior.variance[name]));
    U.$("posterior-mean").textContent = U.number(centre);
    names.forEach((name, i) => { U.$(`${name}-sd`).textContent = `${U.number(sds[i])} km/s`; });
    const limits = [Math.min(centre - 4 * sds[2], 2), Math.max(centre + 4 * sds[2], 2)];
    const density = U.plot(U.$("predictive"), limits, [0.5, 4], "Velocity (km/s)", "");
    density.clip();
    U.line(density, [[2, 0.5], [2, 4]], "#85909F", 1.5, [3, 4]);
    names.forEach((name, i) => {
      const baseline = 3 - i, sd = sds[i];
      const curve = U.linspace(centre - 4 * sd, centre + 4 * sd, 400)
        .map(x => [x, baseline + 0.56 * Math.exp(-0.5 * ((x - centre) / sd) ** 2)]);
      U.area(density, curve, `${colours[i]}33`, baseline);
      U.line(density, curve, colours[i], 2.5);
      U.line(density, posterior.intervals[name].map(x => [x, baseline - 0.09]), colours[i], 4);
      density.ctx.fillStyle = colours[i]; density.ctx.font = "12px system-ui"; density.ctx.textAlign = "left";
      density.ctx.fillText(labels[i], density.left + 5, density.yp(baseline + 0.74));
    });
    density.ctx.restore();
    U.table("intervals", names.map((name, i) => [labels[i], U.number(posterior.variance[name]),
      posterior.intervals[name].map(value => U.number(value)).join(" to ")]));
    if (p.show_data) {
      const observed = U.plot(U.$("observations"), [0, Math.max(2, p.n + 1)], U.extent([...result.x, 2, centre], 0.15),
        "Measured star index", "Velocity (km/s)");
      observed.clip();
      U.line(observed, [[0, 2], [p.n + 1, 2]], "#85909F", 1.5, [3, 4]);
      U.line(observed, [[0, centre], [p.n + 1, centre]], colours[0], 2);
      U.points(observed, result.x.map((x, i) => [i + 1, x]), colours[2], 3);
      observed.ctx.restore();
    } else U.blank(U.$("observations"), "Observations hidden; inference unchanged.");
    U.$("sample-summary").textContent = `n = ${p.n}; observed sample mean ${U.number(result.x.reduce((s, x) => s + x, 0) / p.n)} km/s. ` +
      "Dotted grey: true population mean. Solid purple: posterior mean.";
    const observationVariance = p.intrinsic ** 2 + p.measurement ** 2;
    const sizes = Array.from({ length: 300 }, (_, i) => i + 1);
    const curves = names.map((_, i) => sizes.map(n => [n, Math.sqrt(1 / (1 / 25 + n / observationVariance) +
      (i === 0 ? 0 : i === 1 ? p.intrinsic ** 2 : observationVariance))]));
    const learning = U.plot(U.$("learning"), [1, 300], [0, curves[2][0][1] * 1.08],
      "Number of measured stars", "Standard deviation (km/s)");
    learning.clip();
    curves.forEach((curve, i) => U.line(learning, curve, colours[i], 2.5));
    U.line(learning, [[p.n, 0], [p.n, curves[2][0][1] * 1.08]], "#85909F", 1.5, [3, 4]);
    curves.forEach((_, i) => U.points(learning, [[p.n, sds[i]]], colours[i], 3));
    learning.ctx.restore();
  }
  window.lesson = AstroLesson.mount({
    compute: M.compute, render,
    presets: {
      few: { n: 2 },
      many: { n: 300 },
      "no-scatter": { intrinsic: 0 }
    }
  });
})();
