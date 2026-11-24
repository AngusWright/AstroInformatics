(function (root) {
  "use strict";
  const $ = id => document.getElementById(id);
  function number(value) {
    if (value === null) return "\u2014";
    Astro.check(Number.isFinite(value), "The requested result exceeds numeric precision.");
    return value === 0 ? "0" : Number(value.toPrecision(4)).toString();
  }
  const linspace = (lo, hi, n = 240) => Array.from({ length: n }, (_, i) => lo + (hi - lo) * i / (n - 1));
  function extent(values, pad = 0.06) {
    Astro.check(values.length > 0 && values.every(Number.isFinite), "Plot values must be finite.");
    let lo = Infinity, hi = -Infinity;
    values.forEach(x => { lo = Math.min(lo, x); hi = Math.max(hi, x); });
    const width = hi === lo ? Math.max(1, Math.abs(lo) * 0.1) : hi - lo;
    return hi === lo ? [lo - width / 2, hi + width / 2] : [lo - width * pad, hi + width * pad];
  }
  function ticks(lo, hi) {
    const rough = (hi - lo) / 5, power = 10 ** Math.floor(Math.log10(rough));
    const step = [1, 2, 5, 10].find(n => n * power >= rough) * power;
    const result = [];
    for (let i = Math.ceil(lo / step); i <= Math.floor(hi / step); i++) result.push(i * step);
    return result;
  }
  function tickText(x) {
    if (Math.abs(x) < 1e-14) return "0";
    return Math.abs(x) >= 10000 || Math.abs(x) < 0.01 ? x.toExponential(0) : Number(x.toPrecision(3)).toString();
  }
  function plot(canvas, xRange, yRange, xLabel, yLabel, equal = false, options = {}) {
    const width = Math.max(240, canvas.clientWidth), height = Math.max(180, canvas.clientHeight);
    const ratio = window.devicePixelRatio || 1;
    canvas.width = Math.round(width * ratio); canvas.height = Math.round(height * ratio);
    const ctx = canvas.getContext("2d");
    ctx.setTransform(ratio, 0, 0, ratio, 0, 0);
    ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, width, height);
    const left = 66, right = width - 18, top = 14, bottom = height - 50;
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
    Astro.check([x0, x1, y0, y1, x1 - x0, y1 - y0].every(Number.isFinite) && x1 > x0 && y1 > y0,
      "The requested plot exceeds numeric precision.");
    const xp = x => left + (x - x0) / (x1 - x0) * (right - left);
    const yp = y => bottom - (y - y0) / (y1 - y0) * (bottom - top);
    ctx.font = "12px system-ui, sans-serif"; ctx.fillStyle = "#596a80";
    ctx.strokeStyle = "#9caabd"; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(left, top); ctx.lineTo(left, bottom); ctx.lineTo(right, bottom); ctx.stroke();
    ctx.textAlign = "center";
    (options.xTicks || ticks(x0, x1)).forEach(x => {
      ctx.beginPath(); ctx.moveTo(xp(x), bottom); ctx.lineTo(xp(x), bottom + 5); ctx.stroke();
      ctx.fillText(options.xFormat ? options.xFormat(x) : tickText(x), xp(x), bottom + 21);
    });
    ctx.textAlign = "right";
    ticks(y0, y1).forEach(y => {
      ctx.beginPath(); ctx.moveTo(left - 5, yp(y)); ctx.lineTo(left, yp(y)); ctx.stroke();
      ctx.fillText(tickText(y), left - 9, yp(y) + 4);
    });
    ctx.fillStyle = "#25344b"; ctx.font = "13px system-ui, sans-serif"; ctx.textAlign = "center";
    ctx.fillText(xLabel, (left + right) / 2, height - 8);
    ctx.save(); ctx.translate(15, (top + bottom) / 2); ctx.rotate(-Math.PI / 2);
    ctx.fillText(yLabel, 0, 0); ctx.restore();
    return { ctx, width, height, xp, yp, left, right, top, bottom,
      clip() { ctx.save(); ctx.beginPath(); ctx.rect(left, top, right - left, bottom - top); ctx.clip(); } };
  }
  function line(frame, points, colour, width = 2, dash = []) {
    if (!points.length) return;
    const { ctx, xp, yp } = frame;
    ctx.strokeStyle = colour; ctx.lineWidth = width; ctx.setLineDash(dash); ctx.beginPath();
    points.forEach(([x, y], i) => i ? ctx.lineTo(xp(x), yp(y)) : ctx.moveTo(xp(x), yp(y)));
    ctx.stroke(); ctx.setLineDash([]);
  }
  function area(frame, points, colour, baseline = 0) {
    if (!points.length) return;
    const { ctx, xp, yp } = frame;
    ctx.fillStyle = colour; ctx.beginPath(); ctx.moveTo(xp(points[0][0]), yp(baseline));
    points.forEach(([x, y]) => ctx.lineTo(xp(x), yp(y)));
    ctx.lineTo(xp(points[points.length - 1][0]), yp(baseline)); ctx.closePath(); ctx.fill();
  }
  function points(frame, values, colour, radius = 3) {
    const { ctx, xp, yp } = frame; ctx.fillStyle = colour;
    values.forEach(([x, y]) => { ctx.beginPath(); ctx.arc(xp(x), yp(y), radius, 0, 2 * Math.PI); ctx.fill(); });
  }
  function blank(canvas, message) {
    const width = Math.max(240, canvas.clientWidth), height = Math.max(180, canvas.clientHeight);
    const ratio = window.devicePixelRatio || 1;
    canvas.width = Math.round(width * ratio); canvas.height = Math.round(height * ratio);
    const ctx = canvas.getContext("2d"); ctx.scale(ratio, ratio);
    ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, width, height);
    ctx.fillStyle = "#596a80"; ctx.font = "13px system-ui, sans-serif"; ctx.textAlign = "center";
    ctx.fillText(message, width / 2, height / 2, width - 24);
  }
  function table(id, rows) {
    $(id).replaceChildren(...rows.map(values => {
      const row = document.createElement("tr");
      values.forEach(value => { const cell = document.createElement("td"); cell.textContent = value; row.append(cell); });
      return row;
    }));
  }
  function outputs() {
    document.querySelectorAll("output[for]").forEach(output => { output.value = $(output.getAttribute("for")).value; });
  }
  root.AstroUI = { $, number, linspace, extent, plot, line, area, points, blank, table, outputs };
})(globalThis);
