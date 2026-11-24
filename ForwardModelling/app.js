(function () {
  "use strict";
  const U = AstroUI, M = ForwardModellingModel;
  const percent = value => `${(100 * value).toFixed(0)}%`;
  function cross(frame, point, colour, size) {
    const { ctx, xp, yp } = frame, x = xp(point[0]), y = yp(point[1]);
    ctx.strokeStyle = colour; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(x - size, y); ctx.lineTo(x + size, y);
    ctx.moveTo(x, y - size); ctx.lineTo(x, y + size); ctx.stroke();
  }
  function render(result, parameters, step) {
    const current = M.atFrame(result, step);
    U.$("model-time").textContent = current.time.toFixed(1);
    U.$("member-count").textContent = result.members;
    U.$("outside-median").textContent = percent(current.median);
    U.$("outside-range").textContent = `${percent(current.min)} – ${percent(current.max)}`;
    const frame = U.plot(U.$("encounter-plot"), [-5, 7], [-5, 5], "Host-centred x (model units)", "Host-centred y (model units)", true);
    frame.clip();
    const circle = U.linspace(0, 2 * Math.PI, 161).map(angle => [2.5 * Math.cos(angle), 2.5 * Math.sin(angle)]);
    U.line(frame, circle, "#b7c0cf", 1, [4, 4]);
    current.ensemble.forEach(positions => U.points(frame, positions, "#5059a626", 2));
    U.points(frame, current.nominal, "#087f8c", 3);
    cross(frame, [0, 0], "#25344b", 7);
    cross(frame, current.perturber, "#d56534", 8);
    frame.ctx.restore();
    const spread = U.plot(U.$("spread-plot"), [0, result.steps * result.dt], [0, 1], "Elapsed model time", "Fraction beyond radius 2.5");
    const history = result.summaries.slice(0, step + 1);
    spread.clip();
    spread.ctx.fillStyle = "#5059a626"; spread.ctx.beginPath();
    history.forEach((summary, i) => {
      const x = spread.xp(i * result.dt), y = spread.yp(summary.max);
      if (i) spread.ctx.lineTo(x, y); else spread.ctx.moveTo(x, y);
    });
    for (let i = history.length - 1; i >= 0; i--) spread.ctx.lineTo(spread.xp(i * result.dt), spread.yp(history[i].min));
    spread.ctx.closePath(); spread.ctx.fill();
    U.line(spread, history.map((summary, i) => [i * result.dt, summary.median]), "#5059a6", 2.5);
    U.points(spread, [[current.time, current.median]], "#5059a6", 3.5);
    spread.ctx.restore();
  }
  globalThis.lesson = AstroLesson.mount({
    compute: M.compute, render, steps: M.STEPS,
    presets: {
      deterministic: { impact: 1.5, speed: 1.2, uncertainty: 0 },
      close: { impact: 0.7, speed: 1.2, uncertainty: 0.35 }
    }
  });
})();
