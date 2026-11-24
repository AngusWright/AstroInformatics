(function () {
  "use strict";
  const U = AstroUI, M = MultipleTestingModel;
  const purple = "#5059A6", teal = "#087F8C", orange = "#D56534";
  const pct = value => `${U.number(100 * value)}%`;
  function skyColour(value) {
    if (value === 999) return "#e5e9ef";
    const t = Math.min(1, Math.abs(value) / 3), start = [233, 238, 244];
    const end = value < 0 ? [38, 57, 86] : [248, 185, 145];
    return `rgb(${start.map((channel, i) => Math.round(channel + t * (end[i] - channel))).join(",")})`;
  }
  function render(result, p, frame) {
    const r = M.reveal(result, frame);
    U.$("tested").textContent = `${r.tested} / ${result.cells.length}`;
    U.$("local-p").textContent = U.number(r.probability.single);
    U.$("global-p").textContent = pct(r.probability.global);
    U.$("expected").textContent = U.number(r.probability.expected);
    U.$("discoveries").textContent = `${r.trueDiscoveries} / ${r.falseDiscoveries}`;
    U.$("fdp").textContent = r.falseDiscoveryFraction === null ? "No discoveries" : pct(r.falseDiscoveryFraction);
    U.$("sky-title").textContent = p.signal_fraction === 0 ? "A noise-only sky" : "Sky with optional positive sources";
    const coords = Array.from({ length: p.side }, (_, i) => i + 1);
    const values = coords.map((_, row) => coords.map((__, column) => {
      const i = row * p.side + column;
      return i < r.tested ? result.cells[i].z : 999;
    }));
    const sky = U.plot(U.$("sky"), [0.5, p.side + 0.5], [0.5, p.side + 0.5],
      "Independent sky element", "Row", true);
    sky.clip();
    U.heatmap(sky, coords, coords, values, skyColour);
    const radius = Math.max(2.5, Math.min(10, Math.abs(sky.xp(1) - sky.xp(0)) * 0.43));
    result.cells.slice(0, r.tested).forEach((cell, i) => {
      if (!cell.hit) return;
      sky.ctx.beginPath();
      sky.ctx.arc(sky.xp(i % p.side + 1), sky.yp(Math.floor(i / p.side) + 1), radius, 0, 2 * Math.PI);
      sky.ctx.strokeStyle = cell.truth ? teal : orange; sky.ctx.lineWidth = 2; sky.ctx.stroke();
    });
    sky.ctx.restore();
    const maximumNulls = result.cells.filter(cell => !cell.truth).length;
    const risk = U.plot(U.$("risk"), [0, Math.max(1, maximumNulls)], [0, 1.02],
      "Independent null locations searched", "False-alarm probability");
    risk.clip();
    U.line(risk, U.linspace(0, maximumNulls, 200).map(n =>
      [n, M.searchFalseProbability(p.threshold, Math.round(n), p.tails === "two").global]), purple, 3);
    U.line(risk, [[0, r.probability.single], [maximumNulls, r.probability.single]], orange, 2, [6, 4]);
    U.points(risk, [[r.nullTrials, r.probability.global]], teal, 5);
    risk.ctx.restore();
    U.$("search-summary").textContent =
      `${r.tested} locations searched; ${r.nullTrials} are null. ${r.discoveries} detections: ` +
      `${r.trueDiscoveries} true and ${r.falseDiscoveries} false. ` +
      (r.discoveries ? "The realised false fraction is a one-survey statistic, not the global false-alarm probability."
        : "No discoveries yet: an observed discovery fraction is undefined, even though search risk is well defined.");
  }
  window.lesson = AstroLesson.mount({
    compute: M.compute, render, steps: 100,
    presets: {
      small: { side: 5, threshold: 3, tails: "one", signal_fraction: 0 },
      large: { side: 60, threshold: 3, tails: "one", signal_fraction: 0 },
      sources: { side: 40, threshold: 3, tails: "one", signal_fraction: 0.02, amplitude: 4 }
    }
  });
})();
