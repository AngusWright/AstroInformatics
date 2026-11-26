(function () {
  "use strict";
  const U = AstroUI, M = FourierModel, N = M.SIZE, TAU = 2 * Math.PI;
  const preset = values => ({ ...M.defaults, ...values });
  const presets = {
    hidden: preset({}),
    superposition: preset({ amplitude_a: 1, kx_a: 6, ky_a: 0, phase_a: 0,
      family_b: "wave", amplitude_b: 0.7, noise: 0.6, mode_x: 6, mode_y: 0 }),
    gaussian: preset({ family_a: "gaussian", amplitude_a: 2, width_a: 6, noise: 0,
      remove_mean: false, filter: "all", mode_x: 0, mode_y: 0 }),
    aperture: preset({ family_a: "square", amplitude_a: 1, width_a: 12, noise: 0,
      remove_mean: false, filter: "lowpass", cutoff: 10, mode_x: 0, mode_y: 0 }),
    sources: preset({ family_a: "point", amplitude_a: 4, shift_x_a: -10, family_b: "point",
      amplitude_b: 4, shift_x_b: 10, noise: 0, remove_mean: false, filter: "all", mode_x: 0, mode_y: 0 }),
    noise: preset({ family_a: "none" })
  };
  let cached = null, spectrumFrame = null, current = null;
  const palette = [[19, 32, 55], [66, 61, 125], [28, 146, 153], [250, 203, 95], [255, 247, 218]];
  const mix = (a, b, t) => a.map((value, i) => value + t * (b[i] - value));
  function diverging(value, limit) {
    const fraction = Math.max(-1, Math.min(1, value / limit));
    return mix([248, 249, 250], fraction < 0 ? [38, 88, 135] : [204, 89, 44], Math.abs(fraction));
  }
  function powerColour(value) {
    const position = Math.max(0, Math.min(1, value)) * (palette.length - 1);
    const index = Math.min(palette.length - 2, Math.floor(position));
    return mix(palette[index], palette[index + 1], position - index);
  }
  function phaseColour(angle, brightness) {
    const hue = (angle + Math.PI) / TAU * 6, sector = Math.floor(hue), fraction = hue - sector;
    const low = 0.17, rising = low + (1 - low) * fraction, falling = 1 - (1 - low) * fraction;
    const colours = [[1, rising, low], [falling, 1, low], [low, 1, rising],
      [low, falling, 1], [rising, low, 1], [1, low, falling]];
    return colours[sector % 6].map(value => value * 255 * Math.max(0, Math.min(1, brightness)));
  }
  function image(canvas, colourAt, fourier = false) {
    const width = Math.max(240, canvas.clientWidth), height = Math.max(260, canvas.clientHeight);
    const ratio = devicePixelRatio || 1, side = Math.min(width - 74, height - 64);
    const left = (width - side) / 2 + 12, top = 14, bottom = top + side;
    canvas.width = Math.round(width * ratio); canvas.height = Math.round(height * ratio);
    const ctx = canvas.getContext("2d");
    ctx.setTransform(ratio, 0, 0, ratio, 0, 0); ctx.fillStyle = "#fff"; ctx.fillRect(0, 0, width, height);
    const tile = document.createElement("canvas"); tile.width = N; tile.height = N;
    const buffer = tile.getContext("2d").createImageData(N, N);
    for (let y = 0; y < N; y++) {
      for (let x = 0; x < N; x++) {
        const index = fourier ? M.indexOf(x - N / 2, y - N / 2) : y * N + x;
        const colour = colourAt(index), pixel = ((N - 1 - y) * N + x) * 4;
        buffer.data[pixel] = colour[0]; buffer.data[pixel + 1] = colour[1];
        buffer.data[pixel + 2] = colour[2]; buffer.data[pixel + 3] = 255;
      }
    }
    tile.getContext("2d").putImageData(buffer, 0, 0);
    ctx.imageSmoothingEnabled = false; ctx.drawImage(tile, left, top, side, side);
    const start = fourier ? -N / 2 : 0;
    const xp = value => left + (value - start + 0.5) / N * side;
    const yp = value => bottom - (value - start + 0.5) / N * side;
    const ticks = fourier ? [-64, -32, 0, 32, 63] : [0, 32, 64, 96, 127];
    ctx.font = "11px system-ui, sans-serif"; ctx.fillStyle = "#596a80"; ctx.strokeStyle = "#9caabd"; ctx.lineWidth = 1;
    ctx.strokeRect(left, top, side, side);
    for (const tick of ticks) {
      ctx.beginPath(); ctx.moveTo(xp(tick), bottom); ctx.lineTo(xp(tick), bottom + 4);
      ctx.moveTo(left - 4, yp(tick)); ctx.lineTo(left, yp(tick)); ctx.stroke();
      ctx.textAlign = "center"; ctx.fillText(String(tick), xp(tick), bottom + 18);
      ctx.textAlign = "right"; ctx.fillText(String(tick), left - 8, yp(tick) + 4);
    }
    ctx.fillStyle = "#25344b"; ctx.font = "12px system-ui, sans-serif"; ctx.textAlign = "center";
    ctx.fillText(fourier ? "kx (cycles / image)" : "x (pixels)", left + side / 2, height - 6);
    ctx.save(); ctx.translate(left - 40, top + side / 2); ctx.rotate(-Math.PI / 2);
    ctx.fillText(fourier ? "ky (cycles / image)" : "y (pixels)", 0, 0); ctx.restore();
    return { ctx, left, top, bottom, side, xp, yp };
  }
  function markSpectrum(frame, p) {
    const { ctx, left, top, side, xp, yp } = frame;
    ctx.save(); ctx.beginPath(); ctx.rect(left, top, side, side); ctx.clip();
    if (p.filter === "lowpass" || p.filter === "highpass") {
      ctx.strokeStyle = "#fff"; ctx.lineWidth = 1.5; ctx.setLineDash([5, 4]);
      ctx.beginPath(); ctx.arc(xp(0), yp(0), p.cutoff / N * side, 0, TAU); ctx.stroke(); ctx.setLineDash([]);
    } else {
      const selected = M.indexOf(p.mode_x, p.mode_y), partner = M.conjugateIndex(selected);
      for (const index of new Set([selected, partner])) {
        const x = xp(M.frequency(index % N)), y = yp(M.frequency(Math.floor(index / N)));
        ctx.beginPath(); ctx.arc(x, y, 7, 0, TAU); ctx.strokeStyle = "#132037"; ctx.lineWidth = 4; ctx.stroke();
        ctx.strokeStyle = "#fff"; ctx.lineWidth = 2; ctx.stroke();
      }
    }
    ctx.restore();
  }
  function controls(p) {
    document.querySelectorAll("[data-component]").forEach(node => {
      node.hidden = !node.dataset.families.split(" ").includes(p[`family_${node.dataset.component}`]);
    });
    U.$("cutoff-controls").hidden = !["lowpass", "highpass"].includes(p.filter);
    U.$("spectrum-scale").className = `colour-scale ${p.spectrum === "phase" ? "phase-scale" : "power-scale"}`;
  }
  window.lesson = AstroLesson.mount({
    presets, steps: M.STEPS,
    compute(p) {
      M.validate(p);
      if (!cached || cached.seed !== p.seed) cached = M.prepare(p.seed);
      return M.compute(p, cached);
    },
    render(result, p, step) {
      controls(p); current = M.atFrame(result, step);
      const s = current, number = U.number, ref = result.powerLimit || 1;
      U.$("pixel-ratio").textContent = s.signalToNoise === null ? "No noise" : number(s.signalToNoise);
      U.$("peak-mode").textContent = s.peak ? `(${s.peak.x}, ${s.peak.y})` : "No non-DC power";
      U.$("mode-amplitude").textContent = number(s.amplitude);
      U.$("retained-energy").textContent = s.retainedFraction === null ? "Undefined: zero image" : `${(100 * s.retainedFraction).toFixed(2)}%`;
      U.$("strongest").disabled = s.peak === null;
      image(U.$("image-plot"), i => diverging(s.image[i], result.imageLimit));
      image(U.$("reconstruction-plot"), i => diverging(s.reconstructed[i], result.reconstructionLimit));
      spectrumFrame = image(U.$("spectrum-plot"), i => {
        if (p.spectrum === "phase") return phaseColour(Math.atan2(s.imaginary[i], s.real[i]), Math.sqrt(s.power[i] / ref));
        const value = p.spectrum === "power" ? s.power[i] / ref : Math.log10(1 + 1e5 * s.power[i] / ref) / Math.log10(1 + 1e5);
        return powerColour(value);
      }, true);
      markSpectrum(spectrumFrame, p);
      U.$("image-scale-note").textContent = `Intensity: -${number(result.imageLimit)} to +${number(result.imageLimit)}; zero is white. Current noise SD ${number(s.sigma)}.`;
      U.$("spectrum-note").textContent = p.spectrum === "phase" ?
        "Phase: -180 to +180 degrees; brightness follows amplitude." :
        `${p.spectrum === "logpower" ? "Log-power stretch (five decades)" : "Linear power"}; reference power ${number(result.powerLimit)}.`;
      const label = { pair: `Pair (${p.mode_x}, ${p.mode_y}) and its conjugate`, lowpass: `Low pass: |k| <= ${p.cutoff}`,
        highpass: `High pass: |k| >= ${p.cutoff}`, all: "All Fourier modes" }[p.filter];
      U.$("reconstruction-note").textContent = `${label}. Intensity: -${number(result.reconstructionLimit)} to +${number(result.reconstructionLimit)} (separate contrast).`;
      const row = p.row * N;
      const observed = Array.from(s.image.slice(row, row + N)), truth = Array.from(result.signal.slice(row, row + N));
      const recovered = Array.from(s.reconstructed.slice(row, row + N));
      const profile = U.plot(U.$("profile-plot"), [0, N - 1], U.extent([...observed, ...truth, ...recovered], 0.1),
        "x (pixels)", "Intensity");
      profile.clip();
      U.line(profile, observed.map((value, i) => [i, value]), "#8f9bab", 1);
      U.line(profile, recovered.map((value, i) => [i, value]), "#5059a6", 2);
      U.line(profile, truth.map((value, i) => [i, value]), "#087f8c", 2.5, [6, 4]);
      profile.ctx.restore();
      U.$("profile-note").textContent = `Row y = ${p.row}. Whole-image RMSE versus noiseless input: observed ${number(s.rawRmse)}, reconstructed ${number(s.filteredRmse)}.`;
      const radius = Math.hypot(p.mode_x, p.mode_y);
      const phase = s.amplitude < 1e-12 ? "undefined at zero amplitude" : `${number(s.phase * 180 / Math.PI)} degrees`;
      U.$("mode-note").textContent = `Selected k = (${p.mode_x}, ${p.mode_y}): ${radius === 0 ? "DC (constant offset)" :
        `wavelength ${number(N / radius)} pixels`}; Fourier phase ${phase}. Expected white-noise power per unremoved bin: ${number(s.expectedNoisePower)}. ${p.remove_mean ? "Image mean subtracted." : "Image mean retained."}`;
    }
  });
  U.$("strongest").addEventListener("click", () => {
    if (!lesson.result || !current?.peak) return;
    lesson.set({ mode_x: current.peak.x, mode_y: current.peak.y, filter: "pair" });
  });
  U.$("spectrum-plot").addEventListener("click", event => {
    if (!lesson.result || !spectrumFrame) return;
    const bounds = event.currentTarget.getBoundingClientRect(), x = event.clientX - bounds.left, y = event.clientY - bounds.top;
    const { left, top, bottom, side } = spectrumFrame;
    if (x < left || x >= left + side || y < top || y >= bottom) return;
    lesson.set({ mode_x: Math.floor((x - left) / side * N) - N / 2,
      mode_y: Math.min(N - 1, Math.floor((bottom - y) / side * N)) - N / 2, filter: "pair" });
  });
})();
