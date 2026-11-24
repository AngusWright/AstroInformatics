(function () {
  "use strict";
  const C = window.ConjugateModel, A = window.Astro, U = window.AstroUI, $ = U.$;
  let current = null, parameters = null, role = null, active = "update";
  function sync() {
    const model = $("model").value;
    document.querySelectorAll("[data-models]").forEach(fieldset => {
      fieldset.hidden = !fieldset.dataset.models.split(" ").includes(model);
      fieldset.disabled = fieldset.hidden;
    });
    $("log_n").disabled = $("no_data").checked;
    U.outputs();
  }
  function read() {
    const p = { ...C.DEFAULTS };
    for (const id of Object.keys(p)) {
      const input = $(id);
      if (!input.matches(":disabled")) A.check(input.checkValidity(), `${id}: ${input.validationMessage}`);
      p[id] = id === "model" ? input.value : id === "no_data" ? input.checked : input.valueAsNumber;
    }
    return p;
  }
  function guarded(action) {
    try { action(); }
    catch (error) {
      if (!(error instanceof A.ValidationError)) throw error;
      current = null; role = null;
      $("error").textContent = error.message; $("error").hidden = false;
      for (const id of ["sample-size", "prior-weight", "posterior-mean", "posterior-sd", "information", "data-description"]) $(id).textContent = "";
      $("summary").replaceChildren(); $("sufficiency").replaceChildren();
      document.querySelectorAll("canvas").forEach(canvas => U.blank(canvas, "Correct the settings to continue."));
    }
  }
  function rebuild() {
    sync();
    parameters = read(); current = C.analyse(parameters); role = null;
    $("error").hidden = true;
    render();
  }
  function densityPlot() {
    const distributions = [current.prior, current.likelihood, current.posterior].filter(Boolean);
    const grids = distributions.map(d => U.linspace(0.001, 0.999, 65).map(p => d.quantile(p)));
    const endpoints = grids.flat(), lo = Math.min(...endpoints), hi = Math.max(...endpoints);
    const xs = [...new Set([...U.linspace(lo, hi, 650), ...endpoints,
      ...distributions.flatMap(d => (d.focus || []).flatMap(f => U.linspace(f.location - 4.5 * f.scale, f.location + 4.5 * f.scale, 81)))])]
      .filter(x => x >= lo && x <= hi && (current.model === "normal" || x > 0) && (current.model !== "binomial" || x < 1))
      .sort((a, b) => a - b);
    const curves = [current.prior, current.likelihood, current.posterior].map(d => d ? xs.map(x => [x, d.density(x)]) : []);
    let peak = 0;
    curves.forEach(curve => curve.forEach(([, y]) => { A.check(Number.isFinite(y), "Density cannot be resolved at the selected boundary."); peak = Math.max(peak, y); }));
    const f = U.plot($("density-plot"), [lo, hi], [0, peak * 1.08],
      ["poisson", "gamma"].includes(current.model) ? "Unknown rate theta" : "Unknown parameter theta", "Density");
    f.clip();
    U.area(f, curves[2], "#c7580020");
    U.line(f, curves[0], "#0072b2", 2);
    U.line(f, curves[1], "#596a80", 2, [5, 4]);
    U.line(f, curves[2], "#c75800", 2.5);
    f.ctx.restore();
  }
  function rolePlots() {
    if (!role) role = C.roleCurve(parameters);
    const xr = [0, 4], xs = role.map(r => Math.log10(r.n));
    const yr = U.extent([current.prior.mean, ...role.flatMap(r => [r.lower, r.upper, r.estimate])]);
    const options = { xTicks: [0, 1, 2, 3, 4], xFormat: x => (10 ** x).toLocaleString("en") };
    const f = U.plot($("role-plot"), xr, yr, "Number of observations n", "Unknown parameter theta", false, options);
    f.clip();
    const band = [...role.map((r, i) => [xs[i], r.lower]), ...role.map((r, i) => [xs[i], r.upper]).reverse()];
    f.ctx.fillStyle = "#c7580028"; f.ctx.beginPath();
    band.forEach(([x, y], i) => i ? f.ctx.lineTo(f.xp(x), f.yp(y)) : f.ctx.moveTo(f.xp(x), f.yp(y)));
    f.ctx.closePath(); f.ctx.fill();
    U.line(f, [[0, current.prior.mean], [4, current.prior.mean]], "#0072b2", 1.7, [6, 4]);
    U.line(f, role.map((r, i) => [xs[i], r.estimate]), "#596a80", 1.5, [2, 4]);
    U.line(f, role.map((r, i) => [xs[i], r.mean]), "#c75800", 2.5);
    if (current.n > 0) U.line(f, [[Math.log10(current.n), yr[0]], [Math.log10(current.n), yr[1]]], "#657487", 1, [4, 4]);
    f.ctx.restore();
    const w = U.plot($("weight-plot"), xr, [0, 1], "Number of observations n", "Weight of the prior mean", false, options);
    w.clip();
    U.area(w, role.map((r, i) => [xs[i], r.prior_weight]), "#0072b218");
    U.line(w, role.map((r, i) => [xs[i], r.prior_weight]), "#0072b2", 2.5);
    if (current.n > 0) U.line(w, [[Math.log10(current.n), 0], [Math.log10(current.n), 1]], "#657487", 1, [4, 4]);
    w.ctx.restore();
    $("role-explanation").textContent = (current.n === 0 ? "Current n = 0: posterior equals prior. These curves show hypothetical nonzero sample sizes. " : "") +
      current.weight_description + (current.model === "normal" ?
        " More observations or a smaller known sigma increase data precision; a broader prior decreases prior precision. The posterior approaches Normal(xbar, variance=sigma^2/n)." :
        " This is an exact weighted-mean coefficient, not a probability that the prior is correct. Gamma rate inference depends on the observed sum as well as n.");
  }
  function sufficiency() {
    const p = parameters;
    const r = C.update("normal", 4, p.sample_mean, p.prior_mean, p.prior_sd, { sigma: p.sigma });
    U.table("sufficiency", [[-1, -0.5, 0.5, 1], [-2, -1, 1, 2]].map(offsets =>
      [offsets.map(x => U.number(x + p.sample_mean)).join(", "), 4, U.number(p.sample_mean),
        U.number(r.posterior.mean), U.number(r.posterior.sd)]));
  }
  function render() {
    if (!current) return;
    const r = current, fmt = U.number;
    $("sample-size").textContent = r.n.toLocaleString("en");
    $("prior-weight").textContent = `${fmt(100 * r.prior_weight)}%`;
    $("posterior-mean").textContent = fmt(r.posterior.mean);
    $("posterior-sd").textContent = fmt(r.posterior.sd);
    $("model-description").textContent = C.MODELS[r.model].description;
    $("data-description").textContent = r.data_description;
    $("likelihood-legend").hidden = !r.likelihood;
    $("update-rule").textContent = C.MODELS[r.model].rule;
    $("prior-share").value = r.prior_weight;
    $("weight-description").textContent = r.weight_description;
    $("information").textContent = r.n === 0 ? "No observations: the posterior is exactly the prior. There is no data estimate or likelihood curve." :
      `Posterior mean = ${fmt(100 * r.prior_weight)}% prior mean + ${fmt(100 * (1 - r.prior_weight))}% data estimate (${fmt(r.estimate)}).` +
      (r.model === "normal" ? ` Posterior precision = ${fmt(1 / parameters.prior_sd ** 2)} + ${fmt(r.n / parameters.sigma ** 2)} = ${fmt(1 / r.posterior.sd ** 2)}.` : "");
    U.table("summary", [["Prior", r.prior], ["Posterior", r.posterior]].map(([label, d]) =>
      [label, d.label, fmt(d.mean), fmt(d.sd), `[${fmt(d.quantile(0.025))}, ${fmt(d.quantile(0.975))}]`]));
    if (active === "update") densityPlot();
    if (active === "role") rolePlots();
    if (active === "explain") sufficiency();
  }
  function selectTab(name) {
    active = name;
    for (const id of ["update", "role", "explain"]) {
      $(id).hidden = id !== name;
      $(`tab-${id}`).setAttribute("aria-selected", String(id === name));
      $(`tab-${id}`).tabIndex = id === name ? 0 : -1;
    }
    render();
  }
  const presets = { strong_prior: [0, 0.2], many_data: [4, 0.2], weak_prior: [1, 10], reset: [1, 1] };
  Object.entries(presets).forEach(([id, [logN, sd]]) => $(id).addEventListener("click", () => guarded(() => {
    $("model").value = "normal"; $("no_data").checked = false;
    for (const [key, value] of Object.entries({ log_n: logN, prior_sd: sd, sample_mean: 1, sigma: 2, prior_mean: -1 })) $(key).value = value;
    rebuild();
  })));
  $("controls").addEventListener("submit", event => event.preventDefault());
  $("controls").addEventListener("change", () => guarded(rebuild));
  $("controls").addEventListener("input", U.outputs);
  document.querySelectorAll('[role="tab"]').forEach(button => {
    button.addEventListener("click", () => guarded(() => selectTab(button.getAttribute("aria-controls"))));
    button.addEventListener("keydown", event => {
      const buttons = [...document.querySelectorAll('[role="tab"]')], i = buttons.indexOf(button);
      const destination = { ArrowRight: (i + 1) % 3, ArrowLeft: (i + 2) % 3, Home: 0, End: 2 }[event.key];
      if (destination === undefined) return;
      event.preventDefault(); buttons[destination].focus(); buttons[destination].click();
    });
  });
  window.addEventListener("error", event => {
    $("error").textContent = `Unexpected error: ${event.message}`; $("error").hidden = false;
  });
  new ResizeObserver(() => guarded(render)).observe($("main"));
  if (matchMedia("(max-width: 800px)").matches) $("settings").open = false;
  guarded(rebuild);
  document.documentElement.dataset.ready = "true";
})();
