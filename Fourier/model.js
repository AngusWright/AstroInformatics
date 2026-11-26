(function (root, factory) {
  "use strict";
  const model = factory(typeof module === "object" && module.exports ? require("../TeachingTools/statistics.js") : root.Astro);
  if (typeof module === "object" && module.exports) module.exports = model;
  else root.FourierModel = model;
})(globalThis, function (A) {
  "use strict";
  const SIZE = 128, PIXELS = SIZE * SIZE, STEPS = 40, TAU = 2 * Math.PI;
  const families = ["none", "wave", "gaussian", "square", "disk", "point"];
  const defaults = {
    family_a: "wave", amplitude_a: 0.6, kx_a: 9, ky_a: 5, phase_a: 35, width_a: 8, shift_x_a: 0, shift_y_a: 0,
    family_b: "none", amplitude_b: 0.7, kx_b: -4, ky_b: 10, phase_b: 75, width_b: 8, shift_x_b: 0, shift_y_b: 0,
    noise: 4, seed: 42, remove_mean: true, filter: "pair", mode_x: 9, mode_y: 5,
    cutoff: 12, spectrum: "logpower", row: 64
  };
  const wrap = (value, size = SIZE) => ((value % size) + size) % size;
  const frequency = (index, size = SIZE) => index < size / 2 ? index : index - size;
  const indexOf = (x, y) => wrap(y) * SIZE + wrap(x);
  const conjugateIndex = index => indexOf(-(index % SIZE), -Math.floor(index / SIZE));

  function bounded(value, lo, hi, name, integer = false) {
    A.check(Number.isFinite(value) && value >= lo && value <= hi && (!integer || Number.isInteger(value)),
      `${name} must be ${integer ? "an integer " : ""}between ${lo} and ${hi}.`);
  }
  function validate(p) {
    for (const suffix of ["a", "b"]) {
      A.check(families.includes(p[`family_${suffix}`]), "Choose a valid real-space function.");
      bounded(p[`amplitude_${suffix}`], 0, 5, "Amplitude");
      for (const axis of ["x", "y"]) {
        bounded(p[`k${axis}_${suffix}`], -24, 24, "Wave frequency");
        bounded(p[`shift_${axis}_${suffix}`], -32, 32, "Circular shift", true);
      }
      bounded(p[`phase_${suffix}`], 0, 360, "Phase");
      bounded(p[`width_${suffix}`], 1, 32, "Function size", true);
    }
    bounded(p.seed, 0, 2147483647, "Seed", true);
    bounded(p.noise, 0, 8, "Noise standard deviation");
    bounded(p.mode_x, -SIZE / 2, SIZE / 2 - 1, "Selected x frequency", true);
    bounded(p.mode_y, -SIZE / 2, SIZE / 2 - 1, "Selected y frequency", true);
    bounded(p.cutoff, 0, 91, "Radial cutoff");
    bounded(p.row, 0, SIZE - 1, "Profile row", true);
    A.check(typeof p.remove_mean === "boolean", "Choose whether to subtract the image mean.");
    A.check(["pair", "lowpass", "highpass", "all"].includes(p.filter), "Choose a valid Fourier filter.");
    A.check(["logpower", "power", "phase"].includes(p.spectrum), "Choose a valid spectrum display.");
  }

  function fft1d(real, imaginary, inverse) {
    const n = real.length;
    for (let i = 1, j = 0; i < n; i++) {
      let bit = n >> 1;
      while (j & bit) { j ^= bit; bit >>= 1; }
      j ^= bit;
      if (i < j) {
        [real[i], real[j]] = [real[j], real[i]];
        [imaginary[i], imaginary[j]] = [imaginary[j], imaginary[i]];
      }
    }
    for (let length = 2; length <= n; length *= 2) {
      const angle = (inverse ? 1 : -1) * TAU / length;
      const wr = Math.cos(angle), wi = Math.sin(angle), half = length / 2;
      for (let start = 0; start < n; start += length) {
        let ar = 1, ai = 0;
        for (let offset = 0; offset < half; offset++) {
          const even = start + offset, odd = even + half;
          const br = real[odd] * ar - imaginary[odd] * ai;
          const bi = real[odd] * ai + imaginary[odd] * ar;
          real[odd] = real[even] - br; imaginary[odd] = imaginary[even] - bi;
          real[even] += br; imaginary[even] += bi;
          const next = ar * wr - ai * wi;
          ai = ar * wi + ai * wr; ar = next;
        }
      }
    }
  }
  function transform2d(values, imaginary = null, inverse = false) {
    const valid = array => (Array.isArray(array) || array instanceof Float64Array) && array.every(Number.isFinite);
    A.check(valid(values), "The transform needs a finite real array.");
    const size = Math.sqrt(values.length);
    A.check(Number.isInteger(size) && size >= 2 && size <= SIZE && (size & (size - 1)) === 0,
      `Use a square, power-of-two grid between 2 and ${SIZE} pixels per side.`);
    A.check(imaginary === null || (valid(imaginary) && imaginary.length === values.length), "Real and imaginary grids must match.");
    A.check(typeof inverse === "boolean", "Transform direction must be a boolean.");
    const real = Float64Array.from(values), imag = imaginary === null ? new Float64Array(values.length) : Float64Array.from(imaginary);
    const rowR = new Float64Array(size), rowI = new Float64Array(size);
    for (const columns of [false, true]) {
      for (let row = 0; row < size; row++) {
        for (let column = 0; column < size; column++) {
          const index = columns ? column * size + row : row * size + column;
          rowR[column] = real[index]; rowI[column] = imag[index];
        }
        fft1d(rowR, rowI, inverse);
        for (let column = 0; column < size; column++) {
          const index = columns ? column * size + row : row * size + column;
          real[index] = rowR[column]; imag[index] = rowI[column];
        }
      }
    }
    if (!inverse) {
      for (let i = 0; i < real.length; i++) { real[i] /= values.length; imag[i] /= values.length; }
    }
    return { real, imaginary: imag };
  }

  function component(p, suffix) {
    const image = new Float64Array(PIXELS), family = p[`family_${suffix}`], amplitude = p[`amplitude_${suffix}`];
    if (family === "none" || amplitude === 0) return image;
    const width = p[`width_${suffix}`], phase = p[`phase_${suffix}`] * Math.PI / 180;
    for (let y = 0; y < SIZE; y++) {
      for (let x = 0; x < SIZE; x++) {
        const dx = x - SIZE / 2, dy = y - SIZE / 2, radius2 = dx * dx + dy * dy;
        let value;
        if (family === "wave") value = Math.cos(TAU * (p[`kx_${suffix}`] * x + p[`ky_${suffix}`] * y) / SIZE + phase);
        else if (family === "gaussian") value = Math.exp(-radius2 / (2 * width * width));
        else if (family === "square") value = Math.abs(dx) <= width && Math.abs(dy) <= width ? 1 : 0;
        else if (family === "disk") value = radius2 <= width * width ? 1 : 0;
        else value = dx === 0 && dy === 0 ? 1 : 0;
        image[indexOf(x + p[`shift_x_${suffix}`], y + p[`shift_y_${suffix}`])] = amplitude * value;
      }
    }
    return image;
  }
  function prepare(seed) {
    bounded(seed, 0, 2147483647, "Seed", true);
    const rng = A.createRng(seed);
    return { seed, noise: Float64Array.from({ length: PIXELS }, () => rng.normal()) };
  }
  function makeMask(p) {
    const mask = new Float64Array(PIXELS), selected = indexOf(p.mode_x, p.mode_y), partner = conjugateIndex(selected);
    for (let i = 0; i < PIXELS; i++) {
      const radius = Math.hypot(frequency(i % SIZE), frequency(Math.floor(i / SIZE)));
      mask[i] = p.filter === "all" || (p.filter === "pair" && (i === selected || i === partner)) ||
        (p.filter === "lowpass" && radius <= p.cutoff) || (p.filter === "highpass" && radius >= p.cutoff) ? 1 : 0;
    }
    return mask;
  }
  function filtered(transform, mask) {
    return transform2d(transform.real.map((value, i) => value * mask[i]),
      transform.imaginary.map((value, i) => value * mask[i]), true).real;
  }
  function maximumAbs(values) {
    return values.reduce((maximum, value) => Math.max(maximum, Math.abs(value)), 0);
  }
  function mean(values) {
    const reference = values[0];
    return reference + values.reduce((sum, value) => sum + (value - reference), 0) / values.length;
  }
  function evaluate(result, fraction) {
    const p = result.parameters, sigma = p.noise * fraction;
    const image = new Float64Array(PIXELS), reconstructed = new Float64Array(PIXELS);
    const real = new Float64Array(PIXELS), imaginary = new Float64Array(PIXELS), power = new Float64Array(PIXELS);
    let total = 0, retained = 0, peak = null, peakPower = 0, rawError = 0, filteredError = 0;
    for (let i = 0; i < PIXELS; i++) {
      image[i] = result.signal[i] + sigma * result.noise[i];
      reconstructed[i] = result.filteredSignal[i] + sigma * result.filteredNoise[i];
      real[i] = result.signalTransform.real[i] + sigma * result.noiseTransform.real[i];
      imaginary[i] = result.signalTransform.imaginary[i] + sigma * result.noiseTransform.imaginary[i];
      power[i] = real[i] ** 2 + imaginary[i] ** 2;
      total += power[i]; retained += power[i] * result.mask[i];
      rawError += (image[i] - result.signal[i]) ** 2;
      filteredError += (reconstructed[i] - result.signal[i]) ** 2;
      if (i !== 0 && i <= conjugateIndex(i) && power[i] > peakPower) { peak = i; peakPower = power[i]; }
    }
    const selected = indexOf(p.mode_x, p.mode_y), magnitude = Math.hypot(real[selected], imaginary[selected]);
    return {
      image, reconstructed, real, imaginary, power, sigma,
      peak: peak === null ? null : { x: frequency(peak % SIZE), y: frequency(Math.floor(peak / SIZE)), power: peakPower },
      amplitude: (selected === conjugateIndex(selected) ? 1 : 2) * magnitude,
      phase: magnitude === 0 ? null : Math.atan2(imaginary[selected], real[selected]),
      energy: total, retainedFraction: total === 0 ? null : retained / total,
      rawRmse: Math.sqrt(rawError / PIXELS), filteredRmse: Math.sqrt(filteredError / PIXELS),
      signalToNoise: sigma === 0 ? null : result.signalRms / sigma,
      expectedNoisePower: sigma * sigma / PIXELS
    };
  }
  function compute(parameters, prepared) {
    const p = { ...parameters };
    validate(p);
    const data = prepared || prepare(p.seed);
    A.check(data.seed === p.seed && data.noise instanceof Float64Array && data.noise.length === PIXELS &&
      data.noise.every(Number.isFinite), "Cached noise does not match these settings.");
    const signal = component(p, "a"), second = component(p, "b");
    for (let i = 0; i < PIXELS; i++) signal[i] += second[i];
    const signalMean = mean(signal), noiseMean = mean(data.noise);
    const noise = data.noise.map(value => value - (p.remove_mean ? noiseMean : 0));
    if (p.remove_mean) for (let i = 0; i < PIXELS; i++) signal[i] -= signalMean;
    const signalTransform = transform2d(signal), noiseTransform = transform2d(noise), mask = makeMask(p);
    const result = {
      parameters: p, data, signal, noise, signalMean, noiseMean, signalTransform, noiseTransform, mask,
      filteredSignal: filtered(signalTransform, mask), filteredNoise: filtered(noiseTransform, mask),
      signalRms: Math.sqrt(signal.reduce((sum, value) => sum + value * value, 0) / PIXELS)
    };
    result.full = evaluate(result, 1);
    result.imageLimit = Math.max(maximumAbs(result.full.image), maximumAbs(signal)) || 1;
    result.reconstructionLimit = Math.max(maximumAbs(result.full.reconstructed), maximumAbs(result.filteredSignal)) || 1;
    result.powerLimit = result.full.power.reduce((maximum, value, i) =>
      Math.max(maximum, value, signalTransform.real[i] ** 2 + signalTransform.imaginary[i] ** 2), 0);
    return result;
  }
  function atFrame(result, frame) {
    bounded(frame, 0, STEPS, "Animation frame", true);
    return frame === STEPS ? result.full : evaluate(result, frame / STEPS);
  }
  return { SIZE, PIXELS, STEPS, defaults, families, validate, prepare, transform2d, indexOf, conjugateIndex, frequency, compute, atFrame };
});
