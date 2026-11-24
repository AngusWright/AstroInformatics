(function () {
  "use strict";
  const M = HierarchicalModelsModel, U = AstroUI, A = Astro;
  window.lesson = AstroLesson.mount({
    compute: M.compute, steps: 60,
    presets: { precise: { noise: 0.3, pooling: "partial" }, uncertain: { noise: 2.5, pooling: "partial" } },
    render(result, p, frame) {
      const state = frame === 60 ? result.final : M.atPhase(result, frame / 60), post = result.posterior;
      U.$("shift").textContent = A.quantile(state.estimate.map((x, i) => Math.abs(x - result.observed[i])), 0.5).toFixed(2);
      U.$("population").textContent = post.populationMean === null ? "Not inferred" :
        `${post.populationMean.toFixed(2)} ± ${post.populationSd.toFixed(2)}`;
      U.$("weight").textContent = `${(100 * A.quantile(post.weight, 0.5)).toFixed(1)}%`;
      U.$("assumed-scatter").textContent = p.pooling === "none" ? "Not pooled" :
        p.pooling === "complete" ? "0 (identical stars)" : `${p.scatter.toFixed(1)} km/s`;
      U.$("phase-note").textContent = (frame === 60 ? "Final marginal posterior intervals. " :
        "Illustrative interpolation only: these intermediate intervals are not posterior credible intervals. ") +
        `True population mean: 1.50 km/s. Simulated dispersion: ${p.scatter.toFixed(1)} km/s. ` +
        (p.pooling === "none" ? "No pooling does not infer a shared population mean." :
          "Population and weight metrics, and the table below, always describe the final inference.");
      U.table("stars", result.observed.map((x, i) => [i + 1, x.toFixed(3), result.error[i].toFixed(3),
        post.mean[i].toFixed(3), Math.sqrt(post.variance[i]).toFixed(3), post.weight[i].toFixed(3)]));
      const bounds = [1.5, ...result.truth];
      result.observed.forEach((x, i) => bounds.push(x - 1.96 * result.error[i], x + 1.96 * result.error[i],
        post.mean[i] - 1.96 * Math.sqrt(post.variance[i]), post.mean[i] + 1.96 * Math.sqrt(post.variance[i])));
      const range = U.extent(bounds), f = U.plot(U.$("stars-plot"), range, [0.3, 14.7],
        "Stellar velocity (km/s)", "Star index");
      f.clip();
      U.line(f, [[1.5, 0.3], [1.5, 14.7]], "#087f8c", 2, [6, 4]);
      result.observed.forEach((x, i) => {
        const row = i + 1, estimate = state.estimate[i], sd = state.sd[i];
        U.line(f, [[x - 1.96 * result.error[i], row + 0.15], [x + 1.96 * result.error[i], row + 0.15]], "#b7c0cf", 2);
        U.points(f, [[x, row + 0.15]], "#7f8a9b", 3);
        U.line(f, [[estimate - 1.96 * sd, row - 0.15], [estimate + 1.96 * sd, row - 0.15]], "#5059a6", 2.5);
        U.points(f, [[estimate, row - 0.15]], "#5059a6", 3.5);
        const px = f.xp(result.truth[i]), py = f.yp(row), ctx = f.ctx;
        ctx.fillStyle = "#087f8c"; ctx.beginPath();
        ctx.moveTo(px, py - 4); ctx.lineTo(px + 4, py); ctx.lineTo(px, py + 4); ctx.lineTo(px - 4, py); ctx.closePath(); ctx.fill();
      });
      f.ctx.restore();
    }
  });
})();
