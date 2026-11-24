(function () {
  "use strict";
  const M = SelectionEffectsModel, U = AstroUI;
  const colours = { data: "#d56534", truth: "#087f8c", text: "#25344b" };
  window.lesson = AstroLesson.mount({
    compute: M.compute, steps: 60,
    presets: { shallow: { limit: 18.5 }, deep: { limit: 23 } },
    render(result, p, frame) {
      const state = frame === 60 ? result.final : M.atPhase(result, frame / 60);
      U.$("detected").textContent = `${state.count.toLocaleString()} / 2,400`;
      U.$("missing").textContent = state.missing.toLocaleString();
      U.$("completeness").textContent = `${(100 * state.completeness).toFixed(1)}%`;
      U.$("detected-mean").textContent = state.mean === null ? "No detections" : state.mean.toFixed(2);
      U.$("plot-heading").textContent = `What survives a magnitude limit of ${state.limit.toFixed(1)}?`;
      U.$("phase-note").textContent = `Current limit ${state.limit.toFixed(2)}; target ${p.limit.toFixed(1)}. ` +
        "Playback restricts this same population from magnitude 26; it is not the passage of physical time.";
      const f = U.plot(U.$("selection-plot"), [30, 800], U.extent([...result.luminosity, -2.2, 2.4], 0.03),
        "Distance (Mpc)", "log₁₀(L/L₀)");
      const visible = [], missing = [];
      result.distance.forEach((distance, i) => (state.detected[i] ? visible : missing).push([distance, result.luminosity[i]]));
      f.clip();
      if (p.overlay === "show") U.points(f, missing, "#b7c0cf80", 2);
      U.points(f, visible, "#d5653499", 2.3);
      U.line(f, U.linspace(30, 800).map(x => [x, M.threshold(x, state.limit)]), colours.text, 2);
      U.line(f, [[30, 0], [800, 0]], colours.truth, 2, [6, 4]);
      let segment = [];
      state.binMeans.forEach((mean, i) => {
        if (mean !== null) segment.push([40 + 80 * i, mean]);
        if (mean === null || i === 9) {
          U.line(f, segment, colours.data, 3); U.points(f, segment, colours.data, 3); segment = [];
        }
      });
      f.ctx.restore();
    }
  });
})();
