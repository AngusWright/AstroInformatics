(function (root) {
  "use strict";
  const A = typeof module === "object" && module.exports ? require("../TeachingTools/statistics.js") : root.Astro;
  const rho = 0.45;

  function validate(weight, width) {
    A.check(Number.isFinite(weight) && weight >= 0 && weight <= 1, "Mixture weight must be between 0 and 1.");
    A.check(Number.isFinite(width) && width > 0, "Broad-component SD must be positive.");
  }
  function component(x, y, mx, my, sx, sy) {
    const u = (x - mx) / sx, v = (y - my) / sy;
    return Math.exp(-(u * u - 2 * rho * u * v + v * v) / (2 * (1 - rho * rho))) /
      (2 * Math.PI * sx * sy * Math.sqrt(1 - rho * rho));
  }
  function mixtureJoint(x, y, weight, width) {
    validate(weight, width);
    A.check(Number.isFinite(x) && Number.isFinite(y), "Joint-density coordinates must be finite.");
    return (1 - weight) * component(x, y, -1, 0, 0.23, 0.3) +
      weight * component(x, y, 1.3, 1, width, width);
  }
  function mixtureMarginal(x, weight, width) {
    validate(weight, width);
    A.check(Number.isFinite(x), "Marginal coordinate must be finite.");
    return (1 - weight) * A.normalPdf(x, -1, 0.23) + weight * A.normalPdf(x, 1.3, width);
  }
  function normalInterval(lo, hi, mean, sd) {
    return lo > mean ? A.normalCdf((mean - lo) / sd) - A.normalCdf((mean - hi) / sd) :
      A.normalCdf(hi, mean, sd) - A.normalCdf(lo, mean, sd);
  }
  function mixtureInterval(lower, upper, weight, width) {
    validate(weight, width);
    A.check(Number.isFinite(lower) && Number.isFinite(upper) && upper >= lower, "Interval endpoints must be finite and ordered.");
    return (1 - weight) * normalInterval(lower, upper, -1, 0.23) + weight * normalInterval(lower, upper, 1.3, width);
  }
  function integrateGrid(weight, width, nx = 180, ny = 180) {
    validate(weight, width);
    A.check(Number.isInteger(nx) && nx >= 2 && nx <= 1000 && Number.isInteger(ny) && ny >= 2 && ny <= 1000,
      "Grid dimensions must be integers from 2 to 1000.");
    const dx = 12 / nx, dy = 14 / ny;
    const x = Array.from({ length: nx }, (_, i) => -5 + (i + 0.5) * dx);
    const y = Array.from({ length: ny }, (_, i) => -6 + (i + 0.5) * dy);
    const joint = y.map(v => x.map(u => mixtureJoint(u, v, weight, width)));
    const integrated = x.map((_, i) => joint.reduce((sum, row) => sum + row[i] * dy, 0));
    const mass = integrated.reduce((sum, density) => sum + density * dx, 0);
    return { x, y, joint, integrated, mass, dx, dy };
  }
  function compute(p) {
    A.check(Number.isInteger(p.seed) && p.seed >= 0 && p.seed <= 2147483647, "Seed must be an integer from 0 to 2147483647.");
    A.check(p.weight >= 0.05 && p.weight <= 0.95 && p.width >= 0.5 && p.width <= 2, "Use weight 0.05–0.95 and broad SD 0.5–2.");
    A.check(Number.isFinite(p.centre) && p.centre >= -2 && p.centre <= 3, "Interval centre must be between -2 and 3.");
    A.check(typeof p.show_profile === "boolean", "Profile comparison must be on or off.");
    const grid = integrateGrid(p.weight, p.width);
    const marginal = grid.x.map(x => mixtureMarginal(x, p.weight, p.width));
    const rawProfile = grid.x.map((_, i) => Math.max(...grid.joint.map(row => row[i])));
    const profileArea = rawProfile.reduce((sum, v) => sum + v * grid.dx, 0);
    const lower = p.centre - 0.3, upper = p.centre + 0.3;
    return { ...grid, marginal, profile: rawProfile.map(v => v / profileArea), lower, upper,
      probability: mixtureInterval(lower, upper, p.weight, p.width),
      mode: grid.x[marginal.indexOf(Math.max(...marginal))] };
  }
  const api = { mixtureJoint, mixtureMarginal, mixtureInterval, integrateGrid, compute };
  if (typeof module === "object" && module.exports) module.exports = api;
  else root.MarginalisationModel = api;
})(globalThis);
