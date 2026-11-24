(function () {
  "use strict";
  const M = EddingtonBiasModel, U = AstroUI;
  window.lesson = AstroLesson.mount({
    compute: M.compute, steps: 60,
    presets: { "no-scatter": { noise: 0 }, "steep-noisy": { noise: 0.5, slope: 2.5 } },
    render(result, p, frame) {
      const state = frame === 60 ? result.final : M.atPhase(result, frame / 60), n = result.truth.length;
      U.$("up").textContent = state.up.toLocaleString();
      U.$("down").textContent = state.down.toLocaleString();
      U.$("net").textContent = `${state.net >= 0 ? "+" : ""}${state.net.toLocaleString()}`;
      U.$("observed-high").textContent = state.observedHigh.toLocaleString();
      U.$("phase-note").textContent = `Current error SD ${state.noise.toFixed(3)} dex; target ${p.noise.toFixed(2)} dex. ` +
        `All ${n.toLocaleString()} galaxies are retained, including ${state.observed.filter(x => x < 8.5 || x > 12).length.toLocaleString()} measured outside the true population bounds.`;
      U.table("migration", [
        ["Below threshold", state.lowLow.toLocaleString(), state.up.toLocaleString(), (n - state.trueHigh).toLocaleString()],
        ["At / above threshold", state.down.toLocaleString(), state.highHigh.toLocaleString(), state.trueHigh.toLocaleString()],
        ["Total", (n - state.observedHigh).toLocaleString(), state.observedHigh.toLocaleString(), n.toLocaleString()]
      ]);
      const lo = result.breaks[0], hi = result.breaks[result.breaks.length - 1];
      const grid = [...U.linspace(lo, hi, 650), 8.5 - 1e-8, 8.5, 12, 12 + 1e-8].sort((a, b) => a - b);
      const truth = grid.map(x => [x, M.massDensity(x, p.slope)]);
      const measured = grid.map(x => [x, M.massDensity(x, p.slope, state.noise)]);
      const peak = Math.max(...result.trueHistogram.density, ...state.histogram.density, M.massDensity(8.5, p.slope));
      const f = U.plot(U.$("mass-plot"), [lo, hi], [0, peak * 1.12], "log₁₀ stellar mass (M/M☉)", "Probability density");
      f.clip();
      U.bars(f, result.trueHistogram.mids, result.trueHistogram.density, result.trueHistogram.widths, "#b7c0cf65");
      U.bars(f, state.histogram.mids, state.histogram.density, state.histogram.widths, "#d5653445");
      U.line(f, truth, "#087f8c", 2.5);
      U.line(f, measured, "#d56534", 2.5);
      U.line(f, [[p.threshold, 0], [p.threshold, peak * 1.12]], "#25344b", 2, [6, 4]);
      f.ctx.restore();
    }
  });
})();
