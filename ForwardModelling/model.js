(function (root, factory) {
  "use strict";
  const api = factory(typeof module === "object" && module.exports ? require("../TeachingTools/statistics.js") : root.Astro);
  if (typeof module === "object" && module.exports) module.exports = api;
  root.ForwardModellingModel = api;
})(globalThis, function (A) {
  "use strict";
  const defaults = { seed: 42, impact: 1.5, speed: 1.2, uncertainty: 0.2 };
  const STEPS = 80, DT = 0.1, MEMBERS = 12;

  function validate(p) {
    A.check(Number.isInteger(p.seed) && p.seed >= 0 && p.seed <= 2147483647, "Enter an integer seed between 0 and 2147483647.");
    for (const [key, lo, hi] of [["impact", 0.4, 3], ["speed", 0.6, 3], ["uncertainty", 0, 0.4]]) {
      A.check(Number.isFinite(p[key]) && p[key] >= lo && p[key] <= hi, `${key} must be between ${lo} and ${hi}.`);
    }
  }

  function restrictedEncounter(radii, angles, impact, speed, steps = STEPS, dt = DT) {
    A.check(radii.length > 0 && radii.length === angles.length &&
      radii.every(r => Number.isFinite(r) && r > 0) && angles.every(Number.isFinite), "Provide matching positive radii and finite angles.");
    A.check(Number.isFinite(impact) && impact > 0 && Number.isFinite(speed) && speed > 0, "Encounter parameters must be positive and finite.");
    A.check(Number.isInteger(steps) && steps >= 0 && steps <= 2000 && Number.isFinite(dt) && dt > 0 && dt <= 1,
      "Use 0–2000 integration steps and a time step in (0, 1].");
    let position = radii.map((r, i) => [r * Math.cos(angles[i]), r * Math.sin(angles[i])]);
    const velocity = radii.map((r, i) => {
      const circular = r / (r * r + 0.09) ** 0.75;
      return [-circular * Math.sin(angles[i]), circular * Math.cos(angles[i])];
    });
    const trajectory = [position];
    function acceleration([x, y], time) {
      const px = -5 + speed * time, py = impact;
      const dx = px - x, dy = py - y;
      const host = (x * x + y * y + 0.09) ** 1.5;
      const encounter = (dx * dx + dy * dy + 0.09) ** 1.5;
      const indirect = (px * px + py * py + 0.09) ** 1.5;
      return [-x / host + 0.5 * dx / encounter - 0.5 * px / indirect,
        -y / host + 0.5 * dy / encounter - 0.5 * py / indirect];
    }
    for (let step = 0; step < steps; step++) {
      position.forEach((point, i) => {
        const a = acceleration(point, step * dt);
        velocity[i][0] += dt * a[0] / 2;
        velocity[i][1] += dt * a[1] / 2;
      });
      position = position.map((point, i) => [point[0] + dt * velocity[i][0], point[1] + dt * velocity[i][1]]);
      position.forEach((point, i) => {
        const a = acceleration(point, (step + 1) * dt);
        velocity[i][0] += dt * a[0] / 2;
        velocity[i][1] += dt * a[1] / 2;
      });
      A.check(position.every(point => point.every(Number.isFinite)), "The encounter exceeds numerical precision.");
      trajectory.push(position);
    }
    return trajectory;
  }

  function compute(parameters) {
    const p = { ...defaults, ...parameters };
    validate(p);
    const rng = A.createRng(p.seed);
    const radii = Array.from({ length: 40 }, (_, i) => 0.65 + 0.85 * i / 39);
    const angles = radii.map(() => 2 * Math.PI * rng.uniform());
    const logSd = Math.sqrt(Math.log1p(p.uncertainty ** 2));
    const impacts = Array.from({ length: MEMBERS }, () => p.impact * Math.exp(logSd * rng.normal() - logSd ** 2 / 2));
    const speeds = Array.from({ length: MEMBERS }, () => p.speed * Math.exp(logSd * rng.normal() - logSd ** 2 / 2));
    const nominal = restrictedEncounter(radii, angles, p.impact, p.speed);
    const ensemble = impacts.map((impact, i) => restrictedEncounter(radii, angles, impact, speeds[i]));
    const summaries = nominal.map((_, frame) => {
      const outside = ensemble.map(path => path[frame].filter(([x, y]) => x * x + y * y > 6.25).length / radii.length);
      return { outside, median: A.quantile(outside, 0.5), min: Math.min(...outside), max: Math.max(...outside) };
    });
    return { nominal, ensemble, summaries, impacts, speeds, radii, angles, impact: p.impact, speed: p.speed,
      steps: STEPS, dt: DT, members: MEMBERS };
  }

  function atFrame(result, frame) {
    A.check(Number.isInteger(frame) && frame >= 0 && frame <= result.steps, "Choose an integer encounter frame from 0 to 80.");
    const time = frame * result.dt;
    return { nominal: result.nominal[frame], ensemble: result.ensemble.map(path => path[frame]),
      ...result.summaries[frame], time, perturber: [-5 + result.speed * time, result.impact] };
  }
  return { defaults, STEPS, DT, MEMBERS, validate, restrictedEncounter, compute, atFrame };
});
