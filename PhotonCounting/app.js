(function () {
  "use strict";
  const U = AstroUI, M = PhotonCountingModel;
  let base = null, cached = null, cacheKey = "";
  function stream(canvas, arrivals, elapsed, exposure) {
    const width = Math.max(240, canvas.clientWidth), height = Math.max(180, canvas.clientHeight), ratio = devicePixelRatio || 1;
    canvas.width = Math.round(width * ratio); canvas.height = Math.round(height * ratio);
    const ctx = canvas.getContext("2d"); ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
    ctx.fillStyle = "#ffffff"; ctx.fillRect(0, 0, width, height);
    const left = 40, right = width - 22, bottom = height - 46, top = 22;
    const xp = value => left + value / exposure * (right - left);
    ctx.fillStyle = "#e8f3f3"; ctx.fillRect(left, top, xp(elapsed) - left, bottom - top);
    ctx.strokeStyle = "#9caabd"; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(left, bottom); ctx.lineTo(right, bottom); ctx.stroke();
    ctx.font = "12px system-ui, sans-serif"; ctx.fillStyle = "#596a80"; ctx.textAlign = "center";
    U.linspace(0, exposure, 5).forEach(time => {
      const x = xp(time); ctx.beginPath(); ctx.moveTo(x, bottom); ctx.lineTo(x, bottom + 5); ctx.stroke();
      ctx.fillText(U.number(time), x, bottom + 20);
    });
    ctx.fillStyle = "#25344b"; ctx.font = "13px system-ui, sans-serif";
    ctx.fillText("Time (seconds)", (left + right) / 2, height - 5);
    ctx.strokeStyle = "#d56534"; ctx.lineWidth = 2; ctx.beginPath();
    arrivals.forEach(time => { ctx.moveTo(xp(time), top + 20); ctx.lineTo(xp(time), bottom - 20); }); ctx.stroke();
    ctx.strokeStyle = "#5059a6"; ctx.lineWidth = 2; ctx.setLineDash([6, 4]);
    ctx.beginPath(); ctx.moveTo(xp(elapsed), top); ctx.lineTo(xp(elapsed), bottom); ctx.stroke(); ctx.setLineDash([]);
    if (!arrivals.length) {
      ctx.fillStyle = "#596a80"; ctx.font = "13px system-ui, sans-serif";
      ctx.fillText("No photons received yet", (left + right) / 2, (top + bottom) / 2, right - left - 12);
    }
  }
  function render(result, p, step) {
    const current = M.atFrame(result, step), gaussian = p.gaussian === "show" && current.gaussian !== null;
    U.$("expected-count").textContent = U.number(current.lambda);
    U.$("stream-count").textContent = current.arrivals.length;
    U.$("count-sd").textContent = `${U.number(current.sd)} / ` +
      (current.relativeSd === null ? "undefined" : `${U.number(100 * current.relativeSd)}%`);
    U.$("empirical-moments").textContent = `${U.number(current.mean)} / ${U.number(current.variance)}`;
    U.$("negative-mass").textContent = current.negativeMass === null ? "No Gaussian" : `${U.number(100 * current.negativeMass)}%`;
    U.$("stream-heading").textContent = `One possible photon stream · ${U.number(current.exposure)} of ${U.number(result.exposure)} seconds`;
    stream(U.$("stream-plot"), current.arrivals, current.exposure, result.exposure);
    const maximum = Math.max(...current.empirical, ...current.exact, ...(gaussian ? current.gaussian : []));
    const frame = U.plot(U.$("count-plot"), [current.lower - 0.5, current.upper + 0.5], [0, maximum * 1.12],
      "Photons counted in an exposure", "Probability per integer count");
    frame.clip();
    U.bars(frame, current.empirical.map((_, k) => k), current.empirical, 0.65, "#d5653466");
    if (gaussian) {
      const negative = current.integers.filter(k => k < 0);
      U.bars(frame, negative, negative.map(k => current.gaussian[k - current.lower]), 0.65, "#d5653433");
      U.line(frame, current.integers.map((k, i) => [k, current.gaussian[i]]), "#5059a6", 2.5);
    }
    U.points(frame, current.integers.filter(k => k >= 0).map(k => [k, current.exact[k - current.lower]]), "#087f8c", 3);
    frame.ctx.restore();
    U.$("count-note").textContent = `λ = rate × elapsed exposure = ${U.number(current.lambda)}. ` +
      (current.lambda === 0 ? "Exact point mass at zero; no nondegenerate Gaussian or relative precision is defined. " :
        gaussian ? "The Gaussian is integrated over unit-width integer bins, including impossible negative counts. " :
          "Gaussian curve hidden; its impossible-mass metric remains available. ") +
      "The repeated-exposure variance uses the sample denominator 2,999. Curves retain their original probabilities, without tail renormalisation.";
  }
  globalThis.lesson = AstroLesson.mount({
    compute(parameters) {
      M.validate(parameters);
      if (base === null || base.seed !== parameters.seed) base = M.prepare(parameters.seed);
      const key = `${parameters.seed}/${parameters.rate}/${parameters.exposure}`;
      if (key !== cacheKey) { cached = M.compute(parameters, base); cacheKey = key; }
      return cached;
    },
    render, steps: M.STEPS,
    presets: {
      few: { rate: 0.4, exposure: 1, gaussian: "show" },
      many: { rate: 8, exposure: 15, gaussian: "show" },
      dark: { rate: 0, exposure: 5, gaussian: "show" }
    }
  });
})();
