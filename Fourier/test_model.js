"use strict";
const assert = require("node:assert/strict");
const { test } = require("node:test");
const A = require("../TeachingTools/statistics.js");
const M = require("./model.js");
const base = { ...M.defaults, noise: 0 };
const close = (actual, expected, tolerance = 1e-11) =>
  assert.ok(Math.abs(actual - expected) <= tolerance, `${actual} != ${expected} (tolerance ${tolerance})`);
const maxError = (a, b) => a.reduce((maximum, value, i) => Math.max(maximum, Math.abs(value - b[i])), 0);

test("the complex 2D FFT agrees with a direct DFT, round-trips, and obeys Parseval", () => {
  const n = 8, values = Array.from({ length: n * n }, (_, i) => Math.sin(i * 0.7) + (i % 5) / 3);
  const imaginary = values.map((value, i) => Math.cos(i * 0.3) - value / 2);
  const actual = M.transform2d(values, imaginary);
  for (let ky = 0; ky < n; ky++) for (let kx = 0; kx < n; kx++) {
    let real = 0, imag = 0;
    for (let y = 0; y < n; y++) for (let x = 0; x < n; x++) {
      const angle = -2 * Math.PI * (kx * x + ky * y) / n, i = y * n + x;
      real += values[i] * Math.cos(angle) - imaginary[i] * Math.sin(angle);
      imag += values[i] * Math.sin(angle) + imaginary[i] * Math.cos(angle);
    }
    close(actual.real[ky * n + kx], real / (n * n));
    close(actual.imaginary[ky * n + kx], imag / (n * n));
  }
  const inverse = M.transform2d(actual.real, actual.imaginary, true);
  assert.ok(maxError(inverse.real, values) < 1e-12);
  assert.ok(maxError(inverse.imaginary, imaginary) < 1e-12);
  close(actual.real.reduce((sum, value, i) => sum + value ** 2 + actual.imaginary[i] ** 2, 0),
    values.reduce((sum, value, i) => sum + value ** 2 + imaginary[i] ** 2, 0) / values.length);
});

test("an integer sinusoid has exactly the correct conjugate coefficients, amplitude and phase", () => {
  const r = M.compute({ ...base, amplitude_a: 1.4, phase_a: 30 }), i = M.indexOf(9, 5), j = M.conjugateIndex(i);
  close(r.full.real[i], 0.7 * Math.cos(Math.PI / 6));
  close(r.full.imaginary[i], 0.7 * Math.sin(Math.PI / 6));
  close(r.full.real[j], r.full.real[i]); close(r.full.imaginary[j], -r.full.imaginary[i]);
  close(r.full.amplitude, 1.4); close(r.full.phase, Math.PI / 6);
  close(r.full.retainedFraction, 1);
  assert.ok(maxError(r.full.reconstructed, r.signal) < 1e-12);
  for (let k = 0; k < M.PIXELS; k++) if (k !== i && k !== j) assert.ok(r.full.power[k] < 1e-26);
});

test("superposition adds complex coefficients, including reinforcement and destructive interference", () => {
  const first = M.compute({ ...base, amplitude_a: 1, phase_a: 0 });
  const reinforced = M.compute({ ...base, amplitude_a: 1, phase_a: 0, family_b: "wave",
    amplitude_b: 1, kx_b: 9, ky_b: 5, phase_b: 0 });
  const cancelled = M.compute({ ...reinforced.parameters, phase_b: 180 });
  close(reinforced.full.amplitude, 2);
  close(reinforced.full.energy, 4 * first.full.energy);
  assert.ok(cancelled.full.energy < 1e-26);
  const separate = M.compute({ ...reinforced.parameters, kx_b: -4, ky_b: 10 });
  close(separate.full.energy, 1);
  close(separate.full.retainedFraction, 0.5);
});

test("circular translation preserves power and applies the predicted phase ramp", () => {
  const first = M.compute({ ...base, family_a: "gaussian", width_a: 6, remove_mean: false, filter: "all" });
  const shifted = M.compute({ ...first.parameters, shift_x_a: 7, shift_y_a: -5 });
  assert.ok(maxError(first.full.power, shifted.full.power) < 1e-15);
  const i = M.indexOf(3, 5), angle = -2 * Math.PI * (3 * 7 - 5 * 5) / M.SIZE;
  close(shifted.full.real[i], first.full.real[i] * Math.cos(angle) - first.full.imaginary[i] * Math.sin(angle));
  close(shifted.full.imaginary[i], first.full.real[i] * Math.sin(angle) + first.full.imaginary[i] * Math.cos(angle));
});

test("a point has flat power, and narrowing a Gaussian broadens its Fourier support", () => {
  const point = M.compute({ ...base, family_a: "point", amplitude_a: 3, remove_mean: false, filter: "all" });
  point.full.power.forEach(value => close(value, (3 / M.PIXELS) ** 2, 1e-20));
  close(point.full.energy, 9 / M.PIXELS);
  for (const width of [4, 8]) {
    const r = M.compute({ ...base, family_a: "gaussian", width_a: width, remove_mean: false });
    const ratio = r.full.real[M.indexOf(4, 0)] / r.full.real[0];
    close(ratio, Math.exp(-2 * Math.PI ** 2 * width ** 2 * 4 ** 2 / M.SIZE ** 2), 1e-12);
  }
});

test("all filters are conjugate-symmetric real projections, including DC and Nyquist", () => {
  for (const filter of ["all", "pair", "lowpass", "highpass"]) {
    const r = M.compute({ ...base, family_a: "disk", width_a: 12, noise: 0.4, filter, cutoff: 10 });
    r.mask.forEach((value, i) => assert.equal(value, r.mask[M.conjugateIndex(i)]));
    const reconstructed = M.transform2d(r.full.real.map((value, i) => value * r.mask[i]),
      r.full.imaginary.map((value, i) => value * r.mask[i]), true);
    assert.ok(reconstructed.imaginary.every(value => Math.abs(value) < 1e-12));
    assert.ok(maxError(reconstructed.real, r.full.reconstructed) < 1e-12);
    close(r.full.reconstructed.reduce((sum, value) => sum + value ** 2, 0) / M.PIXELS,
      r.full.energy * r.full.retainedFraction);
    if (filter === "all") assert.ok(maxError(r.full.image, r.full.reconstructed) < 1e-12);
  }
  const dc = M.compute({ ...base, amplitude_a: 2, kx_a: 0, ky_a: 0, phase_a: 0,
    remove_mean: false, mode_x: 0, mode_y: 0 });
  close(dc.full.amplitude, 2); close(dc.full.energy, 4);
  assert.equal(dc.mask.reduce((sum, value) => sum + value, 0), 1);
  const zero = M.compute({ ...dc.parameters, remove_mean: true });
  assert.equal(zero.full.retainedFraction, null); assert.equal(zero.full.peak, null);
  const arbitraryConstant = M.compute({ ...dc.parameters, amplitude_a: 0.6, phase_a: 35, remove_mean: true });
  assert.ok(arbitraryConstant.full.image.every(value => value === 0));
  assert.equal(arbitraryConstant.full.retainedFraction, null);
  const checkerboard = Float64Array.from({ length: M.PIXELS }, (_, i) => i % 2 ? -1 : 1);
  const nyquist = M.compute({ ...base, family_a: "none", noise: 1, mode_x: -64, mode_y: 0 },
    { seed: base.seed, noise: checkerboard });
  close(nyquist.full.amplitude, 1);
  assert.equal(nyquist.mask.reduce((sum, value) => sum + value, 0), 1);
  assert.ok(maxError(nyquist.full.reconstructed, checkerboard) < 1e-12);
});

test("fractional frequencies leak between bins, while a weak coherent wave survives pixel noise", () => {
  const leaked = M.compute({ ...base, kx_a: 9.25 });
  assert.ok(leaked.full.retainedFraction < 0.92);
  const hidden = M.compute(M.defaults);
  assert.deepEqual([hidden.full.peak.x, hidden.full.peak.y], [9, 5]);
  assert.ok(hidden.full.signalToNoise < 0.11);
  assert.ok(hidden.full.rawRmse > 3.8);
  assert.ok(hidden.full.filteredRmse < 0.15);
  close(hidden.full.amplitude, 0.6, 0.15);
  close(hidden.full.expectedNoisePower, 16 / M.PIXELS);
});

test("cached noise, display changes and animation reuse one realisation without mutating it", () => {
  const data = M.prepare(42), original = data.noise.slice(), r = M.compute(M.defaults, data);
  const display = M.compute({ ...M.defaults, spectrum: "phase", row: 21, filter: "all" }, data);
  assert.strictEqual(display.data, data); assert.deepEqual(display.full.image, r.full.image);
  const start = M.atFrame(r, 0), middle = M.atFrame(r, M.STEPS / 2);
  assert.deepEqual(start.image, r.signal);
  assert.ok(middle.image.every((value, i) => Math.abs(value - (start.image[i] + r.full.image[i]) / 2) < 1e-14));
  assert.strictEqual(M.atFrame(r, M.STEPS), r.full);
  assert.deepEqual(data.noise, original);
  assert.deepEqual(M.prepare(42).noise, original);
  assert.notDeepEqual(M.prepare(43).noise, original);
});

test("invalid transform grids, parameters, caches and frames are explicitly rejected", () => {
  for (const patch of [{ family_a: "invalid" }, { noise: -1 }, { amplitude_b: Infinity }, { seed: 0.5 },
    { mode_x: 64 }, { mode_y: -65 }, { filter: "invalid" }, { spectrum: "invalid" }, { row: 128 },
    { width_a: 0 }, { kx_b: 25 }, { shift_x_a: 0.5 }, { cutoff: 92 }, { remove_mean: "yes" }]) {
    assert.throws(() => M.compute({ ...base, ...patch }), A.ValidationError);
  }
  for (const values of [[], [0, 0, 0], Array(9).fill(0), [NaN, 0, 0, 0]]) {
    assert.throws(() => M.transform2d(values), A.ValidationError);
  }
  assert.throws(() => M.transform2d([1, 2, 3, 4], [0]), A.ValidationError);
  assert.throws(() => M.compute(base, M.prepare(0)), A.ValidationError);
  assert.throws(() => M.atFrame(M.compute(base), -1), A.ValidationError);
});
