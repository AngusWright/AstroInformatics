(function () {
  "use strict";
  const M = ShiftModel, U = AstroUI, $ = U.$;
  const colours = ["#187e91", "#d05b32", "#7553ae", "#34834c", "#b63b75", "#826226", "#316bd0", "#8d5342"];
  const dataControls = ["model_type", "noise_model", "noise_min", "noise_max", "n_noise", "num_samples", "sample_min", "sample_max", "seed"];
  const fitControls = ["model_min", "model_max", "nmodels", "weighting", "correct_htrsk"];
  let data = null, groups = null, generation = 0, displayRevision = 0;
  const fmt = value => value === null || !Number.isFinite(value) ? "—" : U.number(value);
  const pair = (a, b, label = "SD") => `${fmt(a)}${b === null ? ` · ${label} unavailable` : ` ± ${fmt(b)}`}`;

  function number(id) {
    const raw = $(id).value;
    Astro.check(raw.trim() !== "" && Number.isFinite(Number(raw)), `Enter a finite value for ${$(id).labels[0].textContent.trim()}.`);
    return Number(raw);
  }

  function readData() {
    return {
      family: $("model_type").value, noise: $("noise_model").value,
      noiseMin: number("noise_min"), noiseMax: number("noise_max"),
      repeats: Math.round(10 ** number("n_noise")), samples: Math.round(10 ** number("num_samples")),
      xmin: number("sample_min"), xmax: number("sample_max"), seed: number("seed")
    };
  }

  function readFit() {
    return { widthMin: number("model_min"), widthMax: number("model_max"), models: number("nmodels"),
      weighting: $("weighting").value, residual: $("correct_htrsk").checked };
  }

  function guard(action) {
    try {
      if (action() === false) return;
      $("error").hidden = true;
      $("error").textContent = "";
    } catch (error) {
      if (!(error instanceof Astro.ValidationError)) throw error;
      $("error").textContent = error.message;
      $("error").hidden = false;
      groups = null;
      U.table("fit-table", []);
      U.blank($("fit-plot"), "Check the settings above to continue.");
      $("status").textContent = "No fit displayed. Correct the invalid settings; no defaults were substituted.";
      $("legend").replaceChildren();
    }
  }

  function updateOutputs() {
    U.outputs();
    $("repeats-label").textContent = `${Math.round(10 ** Number($("n_noise").value)).toLocaleString()} realisations (rounded to an integer).`;
    $("samples-label").textContent = `${Math.round(10 ** Number($("num_samples").value)).toLocaleString()} evenly spaced x positions (rounded).`;
  }

  function generate() {
    const next = M.createData(readData());
    data = next;
    generation++;
    $("realisation").max = String(data.observations.length);
    $("realisation").value = "1";
    $("main").dataset.generation = String(generation);
    $("main").dataset.observations = JSON.stringify([data.observations[0][0], data.observations[0].at(-1), data.observations.at(-1)[0]]);
    $("pending").textContent = "Generated. Fit and display controls keep these exact observations.";
    refit();
  }

  function refit() {
    Astro.check(data !== null, "Generate observations before fitting.");
    const settings = readFit();
    groups = M.compare(data, settings);
    $("main").dataset.fitRevision = String(Number($("main").dataset.fitRevision || 0) + 1);
    const c = data.config;
    $("data-summary").textContent = `${c.samples} x positions × ${c.repeats} realisations · seed ${c.seed}`;
    $("truth-summary").textContent = c.family === "linear" ? "y = x · true intercept 0" : "N(0, 1) · true centroid 0";
    $("weight-summary").textContent = settings.residual ? "Residual heuristic · no analytic intervals"
      : { legacy: "Original 1/SD · not inverse variance", equal: "Equal weights", known: "Known inverse variance · 1/SD²" }[settings.weighting];
    render();
  }

  function selectedIndex() {
    const j = number("realisation") - 1;
    Astro.check(Number.isInteger(j) && j >= 0 && j < data.observations.length, "Choose an existing noise realisation.");
    return j;
  }

  function addLegend(label, colour, dashed = false) {
    const item = document.createElement("span"), swatch = document.createElement("span");
    item.className = "legend-item"; swatch.className = "swatch";
    if (dashed) swatch.style.borderTop = `2px dashed ${colour}`;
    else swatch.style.backgroundColor = colour;
    item.append(swatch, document.createTextNode(label));
    $("legend").append(item);
  }

  function vertical(frame, x, colour, dash = []) {
    U.line(frame, [[x, frame.yMin], [x, frame.yMax]], colour, 1.5, dash);
  }

  function plot(xr, yr, xLabel, yLabel) {
    const frame = U.plot($("fit-plot"), xr, yr, xLabel, yLabel);
    frame.yMin = yr[0]; frame.yMax = yr[1];
    return frame;
  }

  function dataPlot(j) {
    const x = data.x, obs = data.observations[j], sigma = data.sigma;
    const curveX = U.linspace(x[0], x.at(-1), 240);
    const truth = curveX.map(v => data.config.family === "linear" ? v : Astro.normalPdf(v));
    const curves = groups.map(g => {
      const count = Math.min(100, g.fits.length);
      const indices = Array.from({ length: count }, (_, i) => count === 1 ? 0 : Math.round(i * (g.fits.length - 1) / (count - 1)));
      return {
        ensemble: indices.map(i => M.predict(curveX, g.fits[i], g.width, data.config.family)),
        selected: M.predict(curveX, g.fits[j], g.width, data.config.family)
      };
    });
    const allY = [...truth, ...obs.map((y, i) => y - sigma[i]), ...obs.map((y, i) => y + sigma[i])];
    for (const c of curves) { allY.push(...c.selected); for (const ys of c.ensemble) allY.push(...ys); }
    const frame = plot(U.extent(x, 0.02), U.extent(allY), "x", "Observed response");
    frame.clip();
    groups.forEach((g, i) => {
      frame.ctx.globalAlpha = 0.08;
      for (const ys of curves[i].ensemble) U.line(frame, curveX.map((v, k) => [v, ys[k]]), colours[i], 1);
      frame.ctx.globalAlpha = 1;
      U.line(frame, curveX.map((v, k) => [v, curves[i].selected[k]]), colours[i], 2);
    });
    U.line(frame, curveX.map((v, i) => [v, truth[i]]), "#172b3b", 2.5, [6, 4]);
    frame.ctx.strokeStyle = "#697b89";
    frame.ctx.lineWidth = 0.7;
    frame.ctx.beginPath();
    for (let i = 0; i < x.length; i++) {
      const px = frame.xp(x[i]), ya = frame.yp(obs[i] - sigma[i]), yb = frame.yp(obs[i] + sigma[i]);
      frame.ctx.moveTo(px, ya); frame.ctx.lineTo(px, yb);
      frame.ctx.moveTo(px - 2, ya); frame.ctx.lineTo(px + 2, ya);
      frame.ctx.moveTo(px - 2, yb); frame.ctx.lineTo(px + 2, yb);
    }
    frame.ctx.stroke();
    U.points(frame, x.map((v, i) => [v, obs[i]]), "#172b3b", x.length > 500 ? 1.3 : 2.3);
    frame.ctx.restore();
    addLegend("True curve", "#172b3b", true);
    $("plot-note").textContent = `Points and ±1 generating-SD bars: realisation ${j + 1}. Bold coloured curves: its fitted models. Faint curves: ${Math.min(100, data.observations.length)} evenly selected realisations per model (display only; every realisation is fitted).`;
  }

  function biasPlot() {
    const values = groups.flatMap(g => g.fits.map(f => f.shift));
    const [lo, hi] = values.reduce(([a, b], v) => [Math.min(a, v), Math.max(b, v)], [-0.2, 0.2]);
    const bins = 40, step = (hi - lo) / bins;
    const histograms = groups.map(g => {
      const counts = Array(bins).fill(0);
      for (const f of g.fits) counts[Math.min(bins - 1, Math.floor((f.shift - lo) / step))]++;
      return counts.flatMap((n, k) => [[lo + k * step, n / (g.fits.length * step)], [lo + (k + 1) * step, n / (g.fits.length * step)]]);
    });
    const top = Math.max(...histograms.flatMap(p => p.map(v => v[1])));
    const frame = plot([lo, hi], [0, top * 1.12], data.config.family === "linear" ? "Fitted intercept (true 0)" : "Fitted centroid (true 0)", "Histogram density");
    frame.clip();
    vertical(frame, 0, "#172b3b", [5, 4]);
    histograms.forEach((points, i) => U.line(frame, points, colours[i], 2));
    frame.ctx.restore();
    addLegend("True zero shift", "#172b3b", true);
    $("plot-note").textContent = "Common-bin, unit-area histograms of every fitted shift. A single realisation cannot establish a bias or an empirical spread. These are simulation distributions, not parameter posteriors.";
  }

  function trendPlot() {
    const widths = groups.map(g => g.width), values = groups.flatMap(g => g.fits.map(f => f.shift));
    const canTrend = data.observations.length >= 2 && widths[0] !== widths.at(-1);
    const trend = canTrend ? M.biasTrend(groups) : [];
    const yrValues = [...values, 0, ...trend.flatMap(t => [t.mean - M.Z95 * t.sd, t.mean + M.Z95 * t.sd])];
    const frame = plot(U.extent(widths, .12), U.extent(yrValues), data.config.family === "linear" ? "Fixed model slope" : "Fixed model width (SD)", "Fitted shift (true 0)");
    frame.clip();
    U.line(frame, [[widths[0], 0], [widths.at(-1), 0]], "#172b3b", 1, [5, 4]);
    groups.forEach((g, i) => {
      frame.ctx.globalAlpha = .3;
      U.points(frame, g.fits.map(f => [g.width, f.shift]), colours[i], 2.1);
      frame.ctx.globalAlpha = 1;
      U.points(frame, [[g.width, g.summary.mean]], colours[i], 5);
      U.line(frame, [[g.width - .008, g.reference.shift], [g.width + .008, g.reference.shift]], "#172b3b", 3);
    });
    if (canTrend) {
      U.line(frame, trend.map(t => [t.x, t.mean]), "#172b3b", 2);
      for (const sign of [-1, 1]) {
        U.line(frame, trend.map(t => [t.x, t.mean + sign * M.Z95 * t.meanSE]), "#172b3b", 1.3, [5, 4]);
        U.line(frame, trend.map(t => [t.x, t.mean + sign * M.Z95 * t.sd]), "#697b89", 1.3, [2, 4]);
      }
    }
    frame.ctx.restore();
    addLegend("Mean ±1.96 MC SE", "#172b3b", true);
    addLegend("Realisation trend ±1.96 SD", "#697b89", true);
    $("plot-note").textContent = canTrend
      ? "Points: every fitted shift (no random jitter); large points: group means; black ticks: noiseless base-weight references. Solid line: descriptive linear trend. Dashed inner and dotted outer bands respect paired realisations; not classical regression prediction intervals or a guarantee of linearity."
      : "Trend bands need at least two realisations and distinct fixed parameters. Points show individual fits; large points are means and black ticks are noiseless references.";
  }

  function render() {
    updateOutputs();
    if (!data || !groups) return false;
    const j = selectedIndex(), mode = $("plot_type").value, scale = number("annot_cex");
    Astro.check(scale >= .5 && scale <= 2.5, "Annotation scaling must be between 0.5 and 2.5.");
    Astro.check(["data_model", "bias", "correlation"].includes(mode), "Choose a supported figure type.");
    $("legend").replaceChildren();
    const title = { data_model: "Data & Model", bias: "Model Bias", correlation: "Bias–Width / Slope" }[mode];
    $("plot-title").textContent = title;
    $("fit-plot").setAttribute("aria-label", `${title}. Numerical results and uncertainty definitions follow below.`);
    $("plot-title").style.fontSize = `${1 + .2 * scale}rem`;
    $("legend").style.fontSize = `${.65 + .2 * scale}rem`;
    if (mode === "data_model") dataPlot(j);
    else if (mode === "bias") biasPlot();
    else trendPlot();
    groups.forEach((g, i) => addLegend(`${data.config.family === "linear" ? "Slope" : "Width"} ${fmt(g.width)}`, colours[i]));
    U.table("fit-table", groups.map(g => {
      const f = g.fits[j], s = g.summary;
      const coverage = s.intervalCount ? `${(100 * s.coverageTruth).toFixed(1)}% / ${(100 * s.coverageReference).toFixed(1)}% (n=${s.intervalCount})` : "Unavailable (n=0)";
      return [fmt(g.width), pair(f.shift, f.se) + (f.boundary ? " · boundary" : ""),
        f.amplitude === null ? "Not fitted" : pair(f.amplitude, f.amplitudeSE),
        pair(s.mean, s.mcse, "MC SE"), fmt(s.sd), pair(s.reference, s.referenceSE), coverage];
    }));
    const boundaryCount = groups.reduce((s, g) => s + g.summary.boundaryCount, 0);
    const unavailable = groups.reduce((s, g) => s + g.fits.length - g.summary.intervalCount, 0);
    $("status").textContent = `${groups.length} fixed models × ${data.observations.length} paired realisations fitted. Shown: realisation ${j + 1}. ${boundaryCount} boundary fits; ${unavailable} fits without analytic intervals.`;
    $("interval-note").textContent = $("correct_htrsk").checked
      ? "Residual-based correction is a data-dependent heuristic: only empirical scatter and MC SE are reported for these noisy fits. The noiseless reference and its local SD still use the base weights. No heuristic coverage claim is made."
      : `${data.config.family === "linear" ? "Intercept SDs and normal intervals are exact" : "Gaussian centroid/amplitude SDs and normal intervals are local approximations"}, conditional on fixed ${data.config.family === "linear" ? "slope" : "width"}, using the generating noise SDs. Truth is known, not uncertain. Coverage of true zero and of the noiseless reference can differ under misspecification.`;
    $("main").dataset.displayRevision = String(++displayRevision);
  }

  function init() {
    $("controls").addEventListener("submit", event => event.preventDefault());
    dataControls.forEach(id => $(id).addEventListener("change", () => {
      updateOutputs();
      $("pending").textContent = "Data settings changed. Click Generate data to apply; current observations remain fixed.";
    }));
    fitControls.forEach(id => $(id).addEventListener("change", () => guard(refit)));
    ["plot_type", "realisation", "annot_cex"].forEach(id => $(id).addEventListener("change", () => guard(render)));
    $("controls").addEventListener("input", updateOutputs);
    $("generate").addEventListener("click", () => guard(generate));
    new ResizeObserver(() => {
      if (groups && $("error").hidden) guard(render);
    }).observe($("fit-plot"));
    if (window.matchMedia("(max-width: 800px)").matches) $("settings").open = false;
    guard(generate);
    document.documentElement.dataset.ready = "true";
  }
  init();
})();
