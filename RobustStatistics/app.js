(function () {
  "use strict";
  const U = AstroUI, M = RobustStatisticsModel;
  let cached = null, dragFrame = null, dragging = false;
  const presets = { clean: { contaminants: 0 }, extreme: { contaminants: 1, position: 12 }, many: { contaminants: 10, position: 8 } };
  const lesson = AstroLesson.mount({
    presets,
    compute(p) {
      M.validate(p);
      if (!cached || cached.seed !== p.seed) cached = M.prepare(p.seed);
      return M.compute(p, cached);
    },
    render(result, p) {
      const s = result.summaries;
      for (const id of ["mean", "median", "sd", "nmad", "trimmed"]) U.$(id).textContent = U.number(s[id]);
      const limits = U.extent([-12.8, 12.8, ...result.observed, s.mean - s.sd, s.mean + s.sd,
        s.median - s.nmad, s.median + s.nmad], 0.02);
      const plot = U.plot(U.$("sample-plot"), limits, [-0.15, 3.8], "Velocity (km/s)", "");
      plot.ctx.fillStyle = "white"; plot.ctx.fillRect(20, 0, plot.left - 23, plot.bottom + 1);
      plot.clip();
      plot.ctx.fillStyle = "#e8f3f3";
      plot.ctx.fillRect(plot.xp(-1), plot.yp(2.9), plot.xp(1) - plot.xp(-1), plot.yp(2.08) - plot.yp(2.9));
      U.line(plot, [[0, -0.1], [0, 3.5]], "#087f8c", 2, [5, 5]);
      U.points(plot, result.observed.slice(p.contaminants).map((value, i) => [value, result.jitter[i + p.contaminants]]), "#7889a8", 3.8);
      U.points(plot, result.observed.slice(0, p.contaminants).map((value, i) => [value, result.jitter[i]]), "#d56534", 5);
      U.line(plot, [[s.mean - s.sd, 1.55], [s.mean + s.sd, 1.55]], "#d56534", 4);
      U.points(plot, [[s.mean, 1.55]], "#d56534", 6);
      U.line(plot, [[s.median - s.nmad, 0.85], [s.median + s.nmad, 0.85]], "#5059a6", 4);
      U.points(plot, [[s.median, 0.85]], "#5059a6", 6);
      U.points(plot, [[s.trimmed, 0.15]], "#0072b2", 6);
      if (p.contaminants > 0) {
        U.line(plot, [[p.position, 2.78], [p.position, 3.25]], "#d56534", 1.5, [3, 3]);
        U.points(plot, [[p.position, 3.25]], "#d56534", 8);
        plot.ctx.textAlign = "center"; plot.ctx.fillStyle = "#a04a23"; plot.ctx.font = "12px system-ui";
        plot.ctx.fillText("drag", plot.xp(p.position), plot.yp(3.55));
      }
      plot.ctx.restore();
      plot.ctx.textAlign = "left"; plot.ctx.fillStyle = "#596a80"; plot.ctx.font = "12px system-ui";
      for (const [label, y] of [["Observed velocities", 2.96], ["Mean ± SD", 1.86], ["Median ± nMAD", 1.16], ["Trimmed mean", 0.46]]) {
        plot.ctx.fillText(label, plot.left + 5, plot.yp(y));
      }
      dragFrame = { plot, limits, result };
      U.$("sample-plot").style.cursor = p.contaminants > 0 ? "ew-resize" : "default";
      U.$("sample-note").textContent = `${p.contaminants} of ${p.n} measurements replaced (${U.number(100 * p.contaminants / p.n)}%). Trimmed mean removes ${Math.floor(p.n * p.trim)} from each tail. Shading: true central population ±1 SD.`;
      const responses = result.influence;
      const yRange = U.extent([0, ...responses.flatMap(row => [row.mean, row.median, row.trimmed])], 0.15);
      const influence = U.plot(U.$("influence-plot"), [-12, 12], yRange, "Contaminant velocity (km/s)", "Location estimate");
      influence.clip();
      U.line(influence, [[p.position, yRange[0]], [p.position, yRange[1]]], "#9caabd", 1.5, [4, 4]);
      for (const [key, colour] of [["mean", "#d56534"], ["median", "#5059a6"], ["trimmed", "#0072b2"]]) {
        U.line(influence, responses.map(row => [row.position, row[key]]), colour, 2.5);
        U.points(influence, [[p.position, s[key]]], colour, 4);
      }
      influence.ctx.restore();
    }
  });
  window.lesson = lesson;
  const canvas = U.$("sample-plot");
  function local(event) {
    const bounds = canvas.getBoundingClientRect();
    return [event.clientX - bounds.left, event.clientY - bounds.top];
  }
  function move(event) {
    if (!dragging || !lesson.result) return;
    const [x] = local(event), { plot, limits } = dragFrame;
    const raw = limits[0] + (x - plot.left) / (plot.right - plot.left) * (limits[1] - limits[0]);
    const position = Math.round(Math.max(-12, Math.min(12, raw)) * 5) / 5;
    if (Number(U.$("position").value) !== position) {
      U.$("position").value = String(position);
      U.$("position").dispatchEvent(new Event("change", { bubbles: true }));
    }
    event.preventDefault();
  }
  canvas.addEventListener("pointerdown", event => {
    if (!lesson.result || !dragFrame || lesson.parameters.contaminants === 0 || event.button !== 0) return;
    const [x, y] = local(event), { plot, result } = dragFrame, p = lesson.parameters;
    const handles = [[p.position, 3.25], ...result.observed.slice(0, p.contaminants).map((value, i) => [value, result.jitter[i]])];
    if (!handles.some(([value, height]) => Math.hypot(plot.xp(value) - x, plot.yp(height) - y) < 22)) return;
    dragging = true; canvas.setPointerCapture(event.pointerId); move(event);
  });
  canvas.addEventListener("pointermove", move);
  for (const name of ["pointerup", "pointercancel", "lostpointercapture"]) canvas.addEventListener(name, () => { dragging = false; });
})();
