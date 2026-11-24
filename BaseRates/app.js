(function () {
  "use strict";
  const U = AstroUI, M = BaseRatesModel;
  const purple = "#5059A6", teal = "#087F8C", orange = "#D56534";
  const pct = value => value === null ? "Undefined (no flags)" : `${U.number(100 * value)}%`;
  function render(result, p, frame) {
    const r = M.reveal(result, frame), visible = p.reveal === "show", count = r.indices.length;
    U.$("prevalence").textContent = pct(result.prevalence);
    U.$("efficiency").textContent = pct(p.sensitivity);
    U.$("ppv").textContent = pct(result.probability);
    U.$("candidate-count").textContent = `${visible ? r.actual.truePositive : "Hidden"} / ${count}`;
    U.$("realised").textContent = visible ? pct(r.genuineFraction) : "Hidden";
    U.$("scanned").textContent = `${r.scanned} / 10000`;
    const columns = Math.max(5, Math.ceil(Math.sqrt(Math.max(1, count) * 1.6)));
    const rows = Math.max(4, Math.ceil(count / columns));
    const mosaic = U.plot(U.$("mosaic"), [0.5, columns + 0.5], [0.5, rows + 0.5],
      "Packed candidate column", "Packed row", true);
    mosaic.clip();
    const dx = Math.abs(mosaic.xp(1) - mosaic.xp(0)) * 0.72;
    const dy = Math.abs(mosaic.yp(1) - mosaic.yp(0)) * 0.72;
    r.indices.forEach((index, i) => {
      mosaic.ctx.fillStyle = visible ? (result.truth[index] ? teal : orange) : purple;
      mosaic.ctx.fillRect(mosaic.xp(i % columns + 1) - dx / 2,
        mosaic.yp(rows - Math.floor(i / columns)) - dy / 2, dx, dy);
    });
    if (!count) {
      mosaic.ctx.fillStyle = "#596a80"; mosaic.ctx.font = "13px system-ui"; mosaic.ctx.textAlign = "center";
      mosaic.ctx.fillText("No candidates flagged yet", (mosaic.left + mosaic.right) / 2, (mosaic.top + mosaic.bottom) / 2);
    }
    mosaic.ctx.restore();
    U.$("mosaic-legend").textContent = visible ? "Teal: genuine event · Orange: false candidate" : "Purple: flagged; identity unknown";
    const expectedTotal = r.expected.truePositive + r.expected.falsePositive;
    const counts = U.plot(U.$("counts"), [0.4, 2.6], [0, Math.max(1, count, expectedTotal) * 1.14],
      "Scanned survey", "Number of flagged objects", false,
      { xTicks: [1, 2], xFormat: value => value === 1 ? "Expected" : "Realised" });
    counts.clip();
    U.bars(counts, [1], [r.expected.truePositive], 0.62, "#087F8C88");
    U.bars(counts, [1], [expectedTotal], 0.62, "#D5653488", r.expected.truePositive);
    if (visible) {
      U.bars(counts, [2], [r.actual.truePositive], 0.62, teal);
      U.bars(counts, [2], [count], 0.62, orange, r.actual.truePositive);
    } else U.bars(counts, [2], [count], 0.62, purple);
    counts.ctx.restore();
    U.table("outcomes", [
      ["Real event, flagged", U.number(r.expected.truePositive), visible ? r.actual.truePositive : "Hidden"],
      ["Non-event, flagged", U.number(r.expected.falsePositive), visible ? r.actual.falsePositive : "Hidden"],
      ["Real event, missed", U.number(r.expected.falseNegative), visible ? r.actual.falseNegative : "Hidden"],
      ["Non-event, not flagged", U.number(r.expected.trueNegative), visible ? r.actual.trueNegative : "Hidden"]
    ]);
    U.$("bayes-rule").textContent =
      `P(real | flag) = q × s / [q × s + (1 − q) × f] ≈ ${pct(result.probability)}. ` +
      `Expected genuine / all flags after ${r.scanned} objects: approximately ${U.number(r.expected.truePositive)} / ${U.number(expectedTotal)}. ` +
      "Displayed probabilities and counts are rounded for readability; the Bayesian probability uses unrounded inputs.";
  }
  window.lesson = AstroLesson.mount({
    compute: M.compute, render, steps: 100,
    presets: {
      rare: { log_prevalence: -4, sensitivity: 0.95, false_positive: 0.01 },
      common: { log_prevalence: -1, sensitivity: 0.95, false_positive: 0.01 },
      empty: { log_prevalence: -6, sensitivity: 0.95, false_positive: 0 }
    }
  });
})();
