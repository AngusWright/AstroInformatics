(function () {
  "use strict";
  const G = window.Gibbs;
  const $ = id => document.getElementById(id);
  const colours = ["#c75800", "#0072b2"];
  const state = { chain: null, update: 0, playing: false, timer: null, tab: "explore", geometry: null, background: null };
  const line = (a, b, n) => Array.from({ length: n }, (_, i) => a + (b - a) * i / (n - 1));
  const unique = values => [...new Set(values)].sort((a, b) => a - b);
  const number = value => {
    if (value === null) return "\u2014";
    if (!Number.isFinite(value)) throw new G.ValidationError("Settings exceed the numeric plotting range.");
    return value === 0 ? "0" : Number(value.toPrecision(4)).toString();
  };

  function pause() {
    clearTimeout(state.timer);
    state.timer = null;
    state.playing = false;
    $("play").textContent = "Play";
    $("playback-state").textContent = "Paused";
  }

  function showError(error, invalidate = false) {
    pause();
    $("error").textContent = error.message;
    $("error").hidden = false;
    if (invalidate) {
      state.chain = null;
      document.querySelectorAll("canvas").forEach(canvas => blank(canvas, "Correct the settings to continue."));
      $("step-description").textContent = "No valid chain for these settings.";
      $("step-table").replaceChildren();
      $("sample-summary").replaceChildren();
      $("sample-count").textContent = "No valid chain for these settings.";
      $("progress-text").textContent = "Correct the settings to continue.";
      $("progress").value = 0;
    }
    updateButtons();
  }

  function guarded(action, invalidate = false) {
    try {
      action();
    } catch (error) {
      if (!(error instanceof G.ValidationError)) throw error;
      showError(error, invalidate);
    }
  }

  function clearError() { $("error").hidden = true; $("error").textContent = ""; }

  function syncFields() {
    const family = $("family").value;
    document.querySelectorAll("[data-families]").forEach(fieldset => {
      const visible = fieldset.dataset.families.split(" ").includes(family);
      fieldset.hidden = !visible;
      fieldset.disabled = !visible;
    });
    document.querySelectorAll('output[for]').forEach(output => {
      output.value = $(output.getAttribute("for")).value;
    });
    $("mixture-warning").hidden = family !== "mixture";
  }

  function settings() {
    const result = {};
    for (const id of Object.keys(G.DEFAULTS)) {
      const input = $(id);
      if (input.matches(":disabled")) continue;
      if (!input.checkValidity()) {
        throw new G.ValidationError(`${document.querySelector(`label[for="${id}"]`)?.textContent.trim() || id}: ${input.validationMessage}`);
      }
      result[id] = id === "family" ? input.value : input.valueAsNumber;
    }
    return result;
  }

  function makeGeometry(chain) {
    const bounds = [0, 1].map(j => {
      const values = chain.states.map(s => s[j]);
      if (chain.target.limits) values.push(...chain.target.limits[j]);
      const lo = chain.target.positive ? 0 : Math.min(...values);
      const hi = Math.max(...values);
      const pad = (hi - lo) * 0.04;
      const result = [chain.target.positive ? 0 : lo - pad, hi + pad];
      if (!result.every(Number.isFinite) || !(result[1] > result[0])) {
        throw new G.ValidationError("Starting coordinates or parameters are too extreme to plot.");
      }
      return result;
    });
    const axes = bounds.map((b, j) => {
      const sorted = chain.states.map(s => s[j]).sort((a, b) => a - b);
      return unique([
        ...line(...b, 40),
        ...(chain.target.limits ? line(...chain.target.limits[j], 101) : []),
        ...line(0, sorted.length - 1, 41).map(i => sorted[Math.round(i)])
      ]);
    });
    const [xs, ys] = axes;
    const logs = xs.map(x => ys.map(y => chain.target.logDensity(x, y)));
    let peak = -Infinity;
    logs.forEach(row => row.forEach(value => {
      if (Number.isNaN(value) || value === Infinity) throw new G.ValidationError("Density exceeded numeric precision. Try less extreme settings.");
      peak = Math.max(peak, value);
    }));
    if (!Number.isFinite(peak)) throw new G.ValidationError("The target cannot be resolved for these starting coordinates.");
    const density = logs.map(row => row.map(z => Math.exp(z - peak)));
    const segments = [];
    for (const level of [0.01, 0.05, 0.2, 0.5, 0.8]) {
      for (let i = 0; i < xs.length - 1; i++) {
        for (let j = 0; j < ys.length - 1; j++) {
          const vertices = [
            [xs[i], ys[j], density[i][j]], [xs[i + 1], ys[j], density[i + 1][j]],
            [xs[i + 1], ys[j + 1], density[i + 1][j + 1]], [xs[i], ys[j + 1], density[i][j + 1]]
          ];
          const intersections = [];
          for (let edge = 0; edge < 4; edge++) {
            const a = vertices[edge], b = vertices[(edge + 1) % 4];
            if ((a[2] >= level) === (b[2] >= level)) continue;
            const t = (level - a[2]) / (b[2] - a[2]);
            intersections.push([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1])]);
          }
          if (intersections.length === 2) segments.push(intersections);
          if (intersections.length === 4) {
            const centreAbove = vertices.reduce((s, v) => s + v[2], 0) / 4 >= level;
            const pairs = (vertices[0][2] >= level) === centreAbove ? [[0, 1], [2, 3]] : [[3, 0], [1, 2]];
            pairs.forEach(pair => segments.push(pair.map(index => intersections[index])));
          }
        }
      }
    }
    return { bounds, segments };
  }

  function rebuild() {
    pause();
    clearError();
    syncFields();
    const next = G.sampleChain(settings());
    const geometry = makeGeometry(next);
    state.chain = next;
    state.geometry = geometry;
    state.background = null;
    state.update = 0;
    $("burn").max = next.sweeps - 1;
    $("burn").value = Math.min($("burn").valueAsNumber, next.sweeps - 1);
    $("jump_to").max = 2 * next.sweeps;
    $("jump_to").value = 0;
    $("progress").max = 2 * next.sweeps;
    $("target-title").textContent = next.target.label;
    syncFields();
    renderEquations();
    render();
  }

  function updateButtons() {
    for (const id of ["play", "step", "sweep", "rewind", "finish", "jump", "diagnostic-finish"]) {
      $(id).disabled = !state.chain;
    }
    if (state.chain) {
      const done = state.update === 2 * state.chain.sweeps;
      $("step").disabled = done;
      $("sweep").disabled = done;
    }
  }

  function seek(update) {
    if (!state.chain) return;
    if (!Number.isInteger(update) || update < 0 || update > 2 * state.chain.sweeps) {
      throw new G.ValidationError(`Choose an integer update from 0 to ${2 * state.chain.sweeps}.`);
    }
    clearError();
    pause();
    state.update = update;
    render();
  }

  function tick() {
    if (!state.playing || !state.chain) return;
    state.update = Math.min(state.update + 1, 2 * state.chain.sweeps);
    if (state.update === 2 * state.chain.sweeps) pause();
    render();
    if (state.playing) state.timer = setTimeout(() => guarded(tick, true), Number($("delay").value));
  }

  function play() {
    if (!state.chain) return;
    clearError();
    if (state.playing) { pause(); return; }
    if (state.update === 2 * state.chain.sweeps) state.update = 0;
    state.playing = true;
    $("play").textContent = "Pause";
    $("playback-state").textContent = "Playing";
    tick();
  }

  function ticks(lo, hi) {
    const rough = (hi - lo) / 5;
    const power = 10 ** Math.floor(Math.log10(rough));
    const step = [1, 2, 5, 10].find(n => n * power >= rough) * power;
    const result = [];
    for (let i = Math.ceil(lo / step); i <= Math.floor(hi / step); i++) result.push(i * step);
    return result;
  }

  function tickText(x) {
    if (Math.abs(x) < 1e-14) return "0";
    return Math.abs(x) >= 10000 || Math.abs(x) < 0.01 ? x.toExponential(0) : Number(x.toPrecision(3)).toString();
  }

  function plot(canvas, xRange, yRange, xLabel, yLabel, equal = false) {
    const width = Math.max(240, canvas.clientWidth);
    const height = Math.max(180, canvas.clientHeight);
    const ratio = window.devicePixelRatio || 1;
    canvas.width = Math.round(width * ratio);
    canvas.height = Math.round(height * ratio);
    const ctx = canvas.getContext("2d");
    ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
    ctx.fillStyle = "white";
    ctx.fillRect(0, 0, width, height);
    const left = 66, right = width - 18, top = 14, bottom = height - 48;
    let [x0, x1] = xRange, [y0, y1] = yRange;
    if (x0 === x1) { x0 -= 0.5; x1 += 0.5; }
    if (y0 === y1) { y0 -= 0.5; y1 += 0.5; }
    if (equal) {
      const aspect = (right - left) / (bottom - top);
      if ((x1 - x0) / (y1 - y0) < aspect) {
        const centre = (x0 + x1) / 2, half = (y1 - y0) * aspect / 2;
        x0 = centre - half; x1 = centre + half;
      } else {
        const centre = (y0 + y1) / 2, half = (x1 - x0) / aspect / 2;
        y0 = centre - half; y1 = centre + half;
      }
    }
    if (![x0, x1, y0, y1, x1 - x0, y1 - y0].every(Number.isFinite)) {
      throw new G.ValidationError("The requested plot exceeds numeric precision.");
    }
    const xp = x => left + (x - x0) / (x1 - x0) * (right - left);
    const yp = y => bottom - (y - y0) / (y1 - y0) * (bottom - top);
    ctx.font = "12px system-ui, sans-serif";
    ctx.fillStyle = "#596a80";
    ctx.strokeStyle = "#9caabd";
    ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(left, top); ctx.lineTo(left, bottom); ctx.lineTo(right, bottom); ctx.stroke();
    ctx.textAlign = "center";
    ticks(x0, x1).forEach(x => {
      ctx.beginPath(); ctx.moveTo(xp(x), bottom); ctx.lineTo(xp(x), bottom + 5); ctx.stroke();
      ctx.fillText(tickText(x), xp(x), bottom + 21);
    });
    ctx.textAlign = "right";
    ticks(y0, y1).forEach(y => {
      ctx.beginPath(); ctx.moveTo(left - 5, yp(y)); ctx.lineTo(left, yp(y)); ctx.stroke();
      ctx.fillText(tickText(y), left - 9, yp(y) + 4);
    });
    ctx.fillStyle = "#25344b"; ctx.font = "13px system-ui, sans-serif"; ctx.textAlign = "center";
    ctx.fillText(xLabel, (left + right) / 2, height - 7);
    ctx.save(); ctx.translate(15, (top + bottom) / 2); ctx.rotate(-Math.PI / 2);
    ctx.fillText(yLabel, 0, 0); ctx.restore();
    return { ctx, width, height, xp, yp, left, right, top, bottom,
      clip() { ctx.save(); ctx.beginPath(); ctx.rect(left, top, right - left, bottom - top); ctx.clip(); } };
  }

  function blank(canvas, message) {
    const width = Math.max(240, canvas.clientWidth), height = Math.max(180, canvas.clientHeight);
    canvas.width = width * (window.devicePixelRatio || 1);
    canvas.height = height * (window.devicePixelRatio || 1);
    const ctx = canvas.getContext("2d");
    ctx.scale(window.devicePixelRatio || 1, window.devicePixelRatio || 1);
    ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, width, height);
    ctx.fillStyle = "#596a80"; ctx.font = "13px system-ui, sans-serif"; ctx.textAlign = "center";
    ctx.fillText(message, width / 2, height / 2, width - 24);
  }

  function path(frame, points, colour, width = 1.5, dash = []) {
    if (!points.length) return;
    const { ctx, xp, yp } = frame;
    ctx.strokeStyle = colour; ctx.lineWidth = width; ctx.setLineDash(dash); ctx.beginPath();
    points.forEach(([x, y], i) => i === 0 ? ctx.moveTo(xp(x), yp(y)) : ctx.lineTo(xp(x), yp(y)));
    ctx.stroke(); ctx.setLineDash([]);
  }

  function point(frame, value, colour, radius, fill = true) {
    frame.ctx.beginPath();
    frame.ctx.arc(frame.xp(value[0]), frame.yp(value[1]), radius, 0, 2 * Math.PI);
    if (fill) { frame.ctx.fillStyle = colour; frame.ctx.fill(); }
    else { frame.ctx.strokeStyle = colour; frame.ctx.lineWidth = 1.5; frame.ctx.stroke(); }
  }

  function drawJoint(view) {
    const chain = state.chain, canvas = $("joint-plot");
    const f = plot(canvas, ...state.geometry.bounds, ...chain.target.axes, !chain.target.positive);
    const key = `${canvas.width}/${canvas.height}`;
    if (!state.background || state.background.key !== key) {
      f.clip();
      state.geometry.segments.forEach(segment => path(f, segment, "#b8c2d0", 1));
      f.ctx.restore();
      const image = document.createElement("canvas");
      image.width = canvas.width; image.height = canvas.height;
      image.getContext("2d").drawImage(canvas, 0, 0);
      state.background = { key, image };
    } else f.ctx.drawImage(state.background.image, 0, 0, f.width, f.height);
    f.clip();
    const burn = $("burn").valueAsNumber;
    if ($("show_samples").checked) G.samples(chain, state.update).forEach((value, i) =>
      point(f, value, i < burn ? "#9aa2ae70" : "#00826060", 2.2));
    const start = Math.max(0, state.update - $("trail").valueAsNumber);
    path(f, chain.states.slice(start, state.update + 1), "#4d596780", 1.4);
    const initial = chain.states[0], px = f.xp(initial[0]), py = f.yp(initial[1]);
    f.ctx.strokeStyle = "#25344b"; f.ctx.lineWidth = 1.5; f.ctx.beginPath();
    f.ctx.moveTo(px - 5, py - 5); f.ctx.lineTo(px + 5, py + 5);
    f.ctx.moveTo(px - 5, py + 5); f.ctx.lineTo(px + 5, py - 5); f.ctx.stroke();
    const colour = colours[view.coordinate];
    const bounds = state.geometry.bounds;
    const fixed = view.coordinate === 0 ?
      [[bounds[0][0], view.before[1]], [bounds[0][1], view.before[1]]] :
      [[view.before[0], bounds[1][0]], [view.before[0], bounds[1][1]]];
    path(f, fixed, colour, 1, [3, 4]);
    if (view.hasDraw) {
      path(f, [view.before, view.after], colour, 2.5);
      const ax = f.xp(view.after[0]), ay = f.yp(view.after[1]);
      const angle = Math.atan2(ay - f.yp(view.before[1]), ax - f.xp(view.before[0]));
      f.ctx.fillStyle = colour; f.ctx.beginPath(); f.ctx.moveTo(ax, ay);
      f.ctx.lineTo(ax - 8 * Math.cos(angle - 0.4), ay - 8 * Math.sin(angle - 0.4));
      f.ctx.lineTo(ax - 8 * Math.cos(angle + 0.4), ay - 8 * Math.sin(angle + 0.4));
      f.ctx.closePath(); f.ctx.fill();
      point(f, view.before, colour, 5, false);
    }
    point(f, view.after, colour, 5);
    f.ctx.restore();
  }

  function drawConditional(view) {
    const d = view.distribution, j = view.coordinate;
    const lo = Math.min(d.range[0], view.before[j], view.after[j]);
    const hi = Math.max(d.range[1], view.before[j], view.after[j]);
    const xs = unique([...line(lo, hi, 120), ...line(...d.range, 300),
      ...d.focus.flatMap(f => line(f.location - 4 * f.scale, f.location + 4 * f.scale, 81))])
      .filter(x => x >= lo && x <= hi);
    const values = xs.map(x => [x, d.density(x)]);
    const peak = Math.max(...values.map(v => v[1]));
    const f = plot($("conditional-plot"), [lo, hi], [0, peak * 1.08],
      state.chain.target.axes[j], "Conditional density");
    f.clip();
    f.ctx.beginPath(); f.ctx.moveTo(f.xp(lo), f.yp(0));
    values.forEach(([x, y]) => f.ctx.lineTo(f.xp(x), f.yp(y)));
    f.ctx.lineTo(f.xp(hi), f.yp(0)); f.ctx.closePath();
    f.ctx.fillStyle = j === 0 ? "#c7580020" : "#0072b220"; f.ctx.fill();
    path(f, values, colours[j], 2);
    path(f, [[view.before[j], 0], [view.before[j], peak]], "#798798", 1, [4, 4]);
    if (view.hasDraw) path(f, [[view.after[j], 0], [view.after[j], peak]], colours[j], 2);
    f.ctx.restore();
    $("conditional-title").textContent = `${view.hasDraw ? "Last draw" : "Next draw"}: parameter ${j + 1}`;
    $("conditional-description").textContent = `${d.label}. Dashed line: previous value; solid line: new draw.`;
  }

  function tableRows(id, rows) {
    $(id).replaceChildren(...rows.map(values => {
      const row = document.createElement("tr");
      values.forEach(value => {
        const cell = document.createElement("td");
        cell.textContent = value; row.append(cell);
      });
      return row;
    }));
  }

  function drawDiagnostics() {
    const chain = state.chain, burn = $("burn").valueAsNumber;
    const completed = G.samples(chain, state.update), kept = G.samples(chain, state.update, burn);
    $("sample-count").textContent = `${completed.length} completed sweeps; burn-in: ${burn}; retained: ${kept.length}.`;
    tableRows("sample-summary", kept.length < 2 ? [
      ["At least two retained sweeps are needed.", "", "", "", ""]
    ] : [0, 1].map(j => {
      const summary = G.summarize(kept.map(s => s[j]));
      return [chain.target.axes[j], number(summary.mean), number(summary.sd),
        number(chain.target.mean?.[j] ?? null), number(chain.target.sd?.[j] ?? null)];
    }));
    for (let j = 0; j < 2; j++) {
      const trace = [chain.states[0], ...completed].map((s, i) => [i, s[j]]);
      const ys = trace.map(p => p[1]);
      if (chain.target.mean) ys.push(chain.target.mean[j]);
      const extent = Math.max(...ys) - Math.min(...ys);
      const pad = extent === 0 ? 1 : extent * 0.08;
      const f = plot($(`trace-${j}`), [0, Math.max(1, completed.length)],
        [Math.min(...ys) - pad, Math.max(...ys) + pad], "Completed sweep", chain.target.axes[j]);
      f.clip();
      path(f, trace, colours[j], 1.2);
      if (chain.target.mean) path(f, [[0, chain.target.mean[j]], [completed.length, chain.target.mean[j]]], "#657487", 1, [3, 4]);
      path(f, [[burn, Math.min(...ys) - pad], [burn, Math.max(...ys) + pad]], "#657487", 1, [4, 4]);
      f.ctx.restore();
      if (kept.length < 2) {
        blank($(`marginal-${j}`), "Play beyond burn-in to see retained samples.");
        blank($(`acf-${j}`), "At least three retained sweeps are needed.");
        continue;
      }
      const values = kept.map(s => s[j]).sort((a, b) => a - b);
      const lo = values[0], hi = values[values.length - 1];
      if (lo === hi) {
        blank($(`marginal-${j}`), "No variation in retained samples.");
        blank($(`acf-${j}`), "Autocorrelation is undefined without variation.");
        continue;
      }
      const iqr = values[Math.floor(0.75 * (values.length - 1))] - values[Math.floor(0.25 * (values.length - 1))];
      const count = Math.max(5, Math.min(60, iqr > 0 ?
        Math.ceil((hi - lo) / (2 * iqr / Math.cbrt(values.length))) : Math.ceil(Math.sqrt(values.length))));
      const width = (hi - lo) / count;
      const bins = Array(count).fill(0);
      values.forEach(x => { bins[Math.min(count - 1, Math.floor((x - lo) / width))]++; });
      const xlim = chain.target.limits ?
        [Math.min(lo, chain.target.limits[j][0]), Math.max(hi, chain.target.limits[j][1])] : [lo, hi];
      const curve = chain.target.marginal ? line(...xlim, 500).map(x => [x, chain.target.marginal(x, j)]) : [];
      const density = bins.map(n => n / (values.length * width));
      const max = Math.max(...density, ...curve.map(p => p[1]));
      const h = plot($(`marginal-${j}`), xlim, [0, max * 1.08], chain.target.axes[j], "Probability density");
      h.clip(); h.ctx.fillStyle = "#c8d1df";
      density.forEach((value, i) => h.ctx.fillRect(h.xp(lo + i * width), h.yp(value),
        Math.max(0, h.xp(lo + (i + 1) * width) - h.xp(lo + i * width) - 1), h.yp(0) - h.yp(value)));
      path(h, curve, colours[j], 2); h.ctx.restore();
      if (kept.length < 3) { blank($(`acf-${j}`), "At least three retained sweeps are needed."); continue; }
      const correlations = G.acf(kept.map(s => s[j]));
      const a = plot($(`acf-${j}`), [0, Math.max(1, correlations.length - 1)], [-1, 1.05], "Lag (sweeps)", "Autocorrelation");
      a.clip();
      path(a, [[0, 0], [correlations.length - 1, 0]], "#9caabd", 1);
      correlations.forEach((value, lag) => path(a, [[lag, 0], [lag, value]], colours[j], 2));
      a.ctx.restore();
    }
  }

  function renderEquations() {
    const texts = {
      gaussian: ["A single vector observation y, known covariance S, and a flat location prior give the Gaussian posterior shown.",
        "S = [[s1^2, rho*s1*s2], [rho*s1*s2, s2^2]]\ny | theta ~ Gaussian_2(theta, S)\n\nFor coordinate j and the other coordinate k:\nm_j = y_j + rho*s_j/s_k * (theta_k-y_k)\ntheta_j | theta_k, y ~ Normal(mean=m_j, variance=s_j^2*(1-rho^2))"],
      student: ["A single Student-t vector observation with known scale matrix S and a flat location prior. Its covariance is nu/(nu-2) times S; scale is not SD.",
        "y | theta ~ Student-t_2(nu, theta, S)\nm_j = y_j + rho*s_j/s_k * (theta_k-y_k)\ndelta_k = ((theta_k-y_k)/s_k)^2\n\ntheta_j | theta_k, y ~ Student-t(df=nu+1, location=m_j,\n  scale=s_j*sqrt((1-rho^2)*(nu+delta_k)/(nu+1)))"],
      mixture: ["Two Gaussian populations share a covariance matrix S. Centres and weights are specified, not inferred. A full conditional is itself a mixture, with weights determined by the held coordinate.",
        "mu_minus = y - separation*(s1,s2)/2\nmu_plus  = y + separation*(s1,s2)/2\np(theta) = w*N_2(theta; mu_minus,S) + (1-w)*N_2(theta; mu_plus,S)\n\nq = P(population 1 | theta_k)\n  = w*N(theta_k; mu_minus[k],s_k^2) /\n    [w*N(theta_k; mu_minus[k],s_k^2) + (1-w)*N(theta_k; mu_plus[k],s_k^2)]\nm[j,c] = mu_c[j] + rho*s_j/s_k*(theta_k-mu_c[k])\nv_j = s_j^2*(1-rho^2)\n\ntheta_j | theta_k ~ q*N(m[j,minus],v_j) + (1-q)*N(m[j,plus],v_j)\n\nRecompute q before each coordinate draw. At zero separation,\nthe mixture reduces to one Gaussian."],
      poisson: ["A hierarchical posterior for n independent unit-exposure counts with total S. Every Gamma parameterisation uses a rate, not a scale.",
        "x_i | lambda ~ Poisson(lambda)\nlambda | beta ~ Gamma(shape=a, rate=beta)\nbeta ~ Gamma(shape=c, rate=d)\n\nlambda | beta,x ~ Gamma(shape=a+S, rate=n+beta)\nbeta | lambda,x ~ Gamma(shape=a+c, rate=d+lambda)"]
    };
    const [description, equations] = texts[state.chain.target.family];
    const p = document.createElement("p"), pre = document.createElement("pre");
    p.textContent = description; pre.textContent = equations;
    $("model-equations").replaceChildren(p, pre);
  }

  function render() {
    if (!state.chain) return;
    const chain = state.chain, view = G.frame(chain, state.update);
    $("progress").value = state.update;
    $("progress-text").textContent = `Update ${state.update} / ${2 * chain.sweeps} \u00b7 Completed sweeps ${Math.floor(state.update / 2)} / ${chain.sweeps}`;
    const j = view.coordinate, other = 1 - j;
    $("step-description").textContent =
      `${view.hasDraw ? `Last draw: parameter ${j + 1}.` : "Ready: draw parameter 1."} ` +
      `Hold ${chain.target.axes[other]} at ${number(view.before[other])}. ` +
      (view.hasDraw ? `New value: ${number(view.after[j])}. ` : "") +
      (state.update === 2 * chain.sweeps ? "All sweeps complete; Play replays the same chain." :
        `Next: parameter ${state.update % 2 + 1}, conditioned on the current other coordinate.`);
    tableRows("step-table", [0, 1].map(k => [chain.target.axes[k], number(view.before[k]),
      number(view.after[k]), k === j ? (view.hasDraw ? "Sampled" : "Next draw") : "Held fixed"]));
    updateButtons();
    if (state.tab === "explore") { drawJoint(view); drawConditional(view); }
    if (state.tab === "diagnostics") drawDiagnostics();
  }

  function tab(name) {
    state.tab = name;
    for (const id of ["explore", "diagnostics", "explain"]) {
      $(id).hidden = id !== name;
      $(`tab-${id}`).setAttribute("aria-selected", String(id === name));
      $(`tab-${id}`).tabIndex = id === name ? 0 : -1;
    }
    render();
  }

  $("controls").addEventListener("submit", event => event.preventDefault());
  $("controls").addEventListener("change", () => guarded(rebuild, true));
  document.querySelectorAll('input[type="range"]').forEach(input => input.addEventListener("input", syncFields));
  ["burn", "trail", "show_samples"].forEach(id => $(id).addEventListener("change", () => guarded(render)));
  $("delay").addEventListener("change", () => {
    if (state.playing) {
      clearTimeout(state.timer);
      state.timer = setTimeout(() => guarded(tick, true), Number($("delay").value));
    }
  });
  $("play").addEventListener("click", () => guarded(play, true));
  $("step").addEventListener("click", () => guarded(() => seek(state.update + 1)));
  $("sweep").addEventListener("click", () => guarded(() => seek(Math.min(2 * state.chain.sweeps, 2 * (Math.floor(state.update / 2) + 1)))));
  $("rewind").addEventListener("click", () => guarded(() => seek(0)));
  ["finish", "diagnostic-finish"].forEach(id => $(id).addEventListener("click", () => guarded(() => seek(2 * state.chain.sweeps))));
  $("jump").addEventListener("click", () => guarded(() => seek($("jump_to").valueAsNumber)));
  document.querySelectorAll('[role="tab"]').forEach(button => {
    button.addEventListener("click", () => guarded(() => tab(button.getAttribute("aria-controls"))));
    button.addEventListener("keydown", event => {
      const buttons = [...document.querySelectorAll('[role="tab"]')], current = buttons.indexOf(button);
      const destination = { ArrowRight: (current + 1) % 3, ArrowLeft: (current + 2) % 3, Home: 0, End: 2 }[event.key];
      if (destination === undefined) return;
      event.preventDefault(); buttons[destination].focus(); buttons[destination].click();
    });
  });
  document.addEventListener("visibilitychange", () => { if (document.hidden) pause(); });
  window.addEventListener("error", event => showError(new Error(`Unexpected error: ${event.message}`), true));
  new ResizeObserver(() => guarded(render, true)).observe($("main"));
  if (window.matchMedia("(max-width: 800px)").matches) $("settings").open = false;
  guarded(rebuild, true);
  document.documentElement.dataset.ready = "true";
})();
