(function () {
  "use strict";
  const U = AstroUI, M = SurveyShiftModel;
  let base = null;
  const percent = value => `${(100 * value).toFixed(1)}%`;
  function boundary(frame, model, correction, xRange, yRange, colour, dash) {
    const [wx, wy] = model.weight, threshold = model.threshold + wx * correction[0] + wy * correction[1];
    if (Math.abs(wy) > 1e-10 * Math.hypot(wx, wy)) {
      U.line(frame, xRange.map(x => [x, (threshold - wx * x) / wy]), colour, 2, dash);
    } else if (wx !== 0) U.line(frame, yRange.map(y => [threshold / wx, y]), colour, 2, dash);
  }
  function render(result, p) {
    U.$("original-accuracy").textContent = percent(result.originalAccuracy);
    U.$("frozen-accuracy").textContent = percent(result.frozenAccuracy);
    U.$("adapted-accuracy").textContent = percent(result.adaptedAccuracy);
    U.$("sample-count").textContent = `${M.SIZE} / ${M.SIZE}`;
    const targetView = p.view === "target", sample = targetView ? result.target : result.original;
    const prediction = targetView ? result.adaptedPrediction : result.originalPrediction;
    const correction = targetView ? result.correction : [0, 0];
    const all = [...result.training, ...result.original, ...result.target];
    const xRange = U.extent(all.map(row => row[0])), yRange = U.extent(all.map(row => row[1]));
    const frame = U.plot(U.$("classification-plot"), xRange, yRange, "Measured colour 1", "Measured colour 2");
    const xs = U.linspace(...xRange, 65), ys = U.linspace(...yRange, 65);
    const regions = ys.map(y => M.predict(result.model, xs.map(x => [x, y]), correction));
    frame.clip(); U.heatmap(frame, xs, ys, regions, value => value === 0 ? "#e9f4f3" : "#efedf8");
    sample.forEach((row, i) => {
      if (result.labels[i] === 0) U.points(frame, [row], "#087f8c", 3.3);
      else {
        frame.ctx.fillStyle = "#5059a6";
        frame.ctx.fillRect(frame.xp(row[0]) - 3.1, frame.yp(row[1]) - 3.1, 6.2, 6.2);
      }
    });
    frame.ctx.strokeStyle = "#d56534"; frame.ctx.lineWidth = 1.7;
    sample.forEach((row, i) => {
      if (prediction[i] === result.labels[i]) return;
      frame.ctx.beginPath(); frame.ctx.arc(frame.xp(row[0]), frame.yp(row[1]), 5.7, 0, 2 * Math.PI); frame.ctx.stroke();
    });
    if (targetView && p.adaptation !== "none") boundary(frame, result.model, [0, 0], xRange, yRange, "#7f8a9b", [6, 4]);
    boundary(frame, result.model, correction, xRange, yRange, "#25344b", []);
    frame.ctx.restore();
    U.$("plot-heading").textContent = !targetView ? "Independent validation under training conditions" :
      p.adaptation === "none" ? "A frozen classifier in a different survey" : "The same classifier after a colour correction";
    const description = p.adaptation === "none" ? "No target correction." :
      p.adaptation === "known" ? "Known external calibration: subtract the specified simulated offset, not an estimate from target labels." :
        "Unlabelled batch mean alignment: subtract target feature mean minus source training feature mean. This assumes a shared class mixture.";
    U.$("correction-note").textContent = `${description} Applied target correction: ` +
      `(${U.number(result.correction[0])}, ${U.number(result.correction[1])}). ` +
      (!targetView ? "The original-survey plot always uses the uncorrected source rule. " : "") +
      "The source training data and fitted rule remain unchanged.";
    const rows = [
      ["Original survey", result.originalPrediction],
      ["New survey, frozen", result.frozenPrediction],
      ["New survey, selected correction", result.adaptedPrediction]
    ].map(([label, predictions]) => [label, ...[0, 1].map(truth => predictions.filter((value, i) =>
      result.labels[i] === truth && value !== truth).length)]);
    U.table("validation-table", rows);
  }
  globalThis.lesson = AstroLesson.mount({
    compute(parameters) {
      M.validate(parameters);
      if (base === null || base.seed !== parameters.seed) base = M.prepare(parameters.seed);
      return M.compute(parameters, base);
    },
    render,
    presets: {
      matched: { offset: 0, noise: 0, view: "target", adaptation: "none" },
      mismatch: { offset: 1.5, noise: 0.3, view: "target", adaptation: "none" },
      corrected: { offset: 1.5, noise: 0.3, view: "target", adaptation: "known" },
      noisy: { offset: 1.5, noise: 1.2, view: "target", adaptation: "known" }
    }
  });
})();
