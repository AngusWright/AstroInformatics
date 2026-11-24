"use strict";
const assert = require("node:assert/strict");
const A = require("../TeachingTools/statistics.js"), M = require("./model.js");
const close = (actual, expected, tolerance = 1e-11) => assert.ok(Math.abs(actual - expected) <= tolerance,
  `${actual} != ${expected} (tolerance ${tolerance})`);
const originalRandom = Math.random;
Math.random = () => { throw new Error("Use the seeded RNG, not Math.random"); };
try {
  // Calibrated against retained model.R: restricted_encounter(c(.65,1,1.5), c(0,pi/3,pi),1.5,1.2).
  const fixture = M.restrictedEncounter([0.65, 1, 1.5], [0, Math.PI / 3, Math.PI], 1.5, 1.2);
  const first = [[0.641158869729952, 0.107303148040639], [0.416627733101116, 0.90907397439834],
    [-1.49797431092862, -0.0792418130140672]];
  const final = [[0.593080275700244, 0.320613945907223], [-0.420207334282248, 0.902198743782879],
    [1.05091677007225, 1.11842473064881]];
  first.forEach((point, i) => point.forEach((value, j) => close(fixture[1][i][j], value)));
  final.forEach((point, i) => point.forEach((value, j) => close(fixture[80][i][j], value)));

  const deterministic = M.compute({ impact: 2, speed: 1, uncertainty: 0 });
  deterministic.ensemble.forEach(path => assert.deepEqual(path, deterministic.nominal));
  assert.ok(deterministic.impacts.every(value => value === 2));
  assert.ok(deterministic.speeds.every(value => value === 1));
  const uncertain = M.compute({ impact: 0.7, uncertainty: 0.35 });
  // R qlnorm(pnorm(z), log(mean)-log1p(.35^2)/2, sqrt(log1p(.35^2))) for the fixed seeded z.
  [0.466805379354913, 0.420478217872444, 0.541632497263592].forEach((value, i) => close(uncertain.impacts[i], value));
  [1.02040322072531, 0.85666025663326, 1.26683111681597].forEach((value, i) => close(uncertain.speeds[i], value));
  assert.deepEqual(uncertain, M.compute({ impact: 0.7, uncertainty: 0.35 }), "Seeded simulations must be repeatable");
  assert.notDeepEqual(uncertain.angles, M.compute({ impact: 0.7, uncertainty: 0.35, seed: 43 }).angles);
  assert.deepEqual(uncertain.nominal, M.compute({ impact: 0.7, uncertainty: 0 }).nominal,
    "Uncertainty must not move the nominal prediction");
  uncertain.ensemble.forEach(path => {
    assert.deepEqual(path[0], uncertain.nominal[0], "Every member must share the same initial disc");
    assert.ok(path.flat(2).every(Number.isFinite));
  });
  assert.ok(uncertain.impacts.every(value => value > 0) && uncertain.speeds.every(value => value > 0));
  assert.ok(new Set(uncertain.impacts).size > 1 && new Set(uncertain.speeds).size > 1);
  const initial = M.atFrame(uncertain, 0), finalFrame = M.atFrame(uncertain, 80);
  assert.equal(initial.time, 0); assert.ok(initial.outside.every(value => value === 0));
  assert.equal(finalFrame.time, 8); close(finalFrame.perturber[0], 4.6); close(finalFrame.perturber[1], 0.7);
  assert.ok(finalFrame.max > finalFrame.min, "Input uncertainty must propagate into the morphological diagnostic");
  assert.equal(finalFrame.nominal, uncertain.nominal[80], "Animation must reveal cached trajectories");
  assert.equal(M.atFrame(uncertain, 10).ensemble[0], uncertain.ensemble[0][10]);
  close(finalFrame.median, 0.2375, 1e-14);
  uncertain.summaries.forEach(summary => {
    assert.ok(summary.min >= 0 && summary.max <= 1 && summary.median >= summary.min && summary.median <= summary.max);
  });

  // A distant perturber approaches the analytic softened circular orbit; halve dt at equal physical time.
  const coarse = M.restrictedEncounter([0.8, 1.2], [0, 1], 1000, 1, 80, 0.05);
  const fine = M.restrictedEncounter([0.8, 1.2], [0, 1], 1000, 1, 160, 0.025);
  for (let i = 0; i < 2; i++) {
    fine[160][i].forEach((value, j) => close(coarse[80][i][j], value, 0.01));
    close(Math.hypot(...fine[160][i]), [0.8, 1.2][i], 0.002);
  }
  for (const impact of [0.4, 3]) for (const speed of [0.6, 3]) for (const uncertainty of [0, 0.4]) {
    const result = M.compute({ seed: impact === 0.4 ? 0 : 2147483647, impact, speed, uncertainty });
    assert.ok(result.nominal.flat(2).every(Number.isFinite));
    assert.ok(result.ensemble.every(path => path.flat(2).every(Number.isFinite)));
  }
  for (const bad of [{ seed: -1 }, { seed: 1.5 }, { seed: 2147483648 }, { impact: 0 }, { impact: NaN },
    { speed: Infinity }, { speed: 0.5 }, { uncertainty: -0.1 }, { uncertainty: 0.41 }]) {
    assert.throws(() => M.compute(bad), A.ValidationError);
  }
  assert.throws(() => M.atFrame(uncertain, 81), A.ValidationError);
  assert.throws(() => M.atFrame(uncertain, 0.5), A.ValidationError);
  assert.throws(() => M.restrictedEncounter([0], [0], 1, 1), A.ValidationError);
  assert.throws(() => M.restrictedEncounter([1], [], 1, 1), A.ValidationError);
  assert.throws(() => M.restrictedEncounter([1], [0], 1, 1, 80, 0), A.ValidationError);
  console.log("ForwardModelling: retained-R trajectories, uncertainty propagation, circular-orbit limits and validation passed.");
} finally { Math.random = originalRandom; }
