"use strict";
const assert = require("node:assert/strict");
const { runLessonBrowser } = require("../TeachingTools/lesson_browser_test.js");

runLessonBrowser("Fourier", async b => {
  assert.equal(await b.text("peak-mode"), "(9, 5)");
  assert.ok(await b.evaluate("lesson.result.full.signalToNoise < 0.11 && lesson.result.full.filteredRmse < 0.15"));
  assert.ok(await b.hasColour("profile-plot", [8, 127, 140]));
  assert.ok(await b.hasColour("profile-plot", [80, 89, 166]));
  await b.change("noise", 0);
  assert.ok(await b.evaluate("lesson.result.full.image.every((x,i) => Math.abs(x - lesson.result.full.reconstructed[i]) < 1e-12)"));
  assert.ok(await b.hasColour("profile-plot", [8, 127, 140]));
  assert.ok(await b.hasColour("profile-plot", [80, 89, 166]));
  await b.click("reset");
  await b.evaluate("window.originalNoise = lesson.result.data; window.originalImage = lesson.result.full.image.slice()");
  for (const display of ["power", "phase", "logpower"]) {
    await b.change("spectrum", display); await b.healthy();
    assert.ok(await b.evaluate("lesson.result.data === originalNoise && lesson.result.full.image.every((x,i) => x === originalImage[i])"));
  }
  await b.change("filter", "all");
  assert.ok(await b.evaluate("lesson.result.full.image.every((x,i) => Math.abs(x - lesson.result.full.reconstructed[i]) < 1e-12)"));
  await b.change("noise", 0);
  for (const family of ["wave", "gaussian", "square", "disk", "point", "none"]) {
    await b.change("family_a", family); await b.healthy();
    assert.ok(await b.evaluate("lesson.result.full.power.every(Number.isFinite)"));
    assert.equal(await b.evaluate('document.querySelector(\'[data-component="a"][data-families="wave"]\').hidden'), family !== "wave");
  }
  assert.equal(await b.text("retained-energy"), "Undefined: zero image");
  await b.change("mode_x", 64);
  assert.equal(await b.evaluate('document.getElementById("error").hidden'), false);
  assert.ok(await b.evaluate('[...document.querySelectorAll("[data-result]")].every(node => node.textContent === "\\u2014")'));
  await b.click("reset"); await b.healthy();

  await b.evaluate('document.querySelector(\'[data-preset="superposition"]\').click()');
  const spot = await b.evaluate(`(() => {
    const canvas = document.getElementById("spectrum-plot");
    canvas.scrollIntoView({ block: "center" });
    const box = canvas.getBoundingClientRect();
    const side = Math.min(canvas.clientWidth - 74, canvas.clientHeight - 64);
    const left = (canvas.clientWidth - side) / 2 + 12;
    return { x: box.left + left + (-4 + 64 + 0.5) / 128 * side,
      y: box.top + 14 + side - (10 + 64 + 0.5) / 128 * side };
  })()`);
  await b.call("Input.dispatchMouseEvent", { type: "mousePressed", ...spot, button: "left", clickCount: 1 });
  await b.call("Input.dispatchMouseEvent", { type: "mouseReleased", ...spot, button: "left", clickCount: 1 });
  assert.equal(await b.value("mode_x"), "-4"); assert.equal(await b.value("mode_y"), "10");
  assert.equal(await b.value("filter"), "pair");
  assert.ok(await b.evaluate("Math.abs(lesson.result.full.amplitude - 0.7) < 0.1"));
  await b.change("filter", "highpass"); await b.change("cutoff", 91); await b.healthy();
  assert.ok(await b.evaluate("lesson.result.full.reconstructed.every(x => x === 0)"));
  await b.change("filter", "lowpass");
  assert.equal(await b.evaluate('document.getElementById("cutoff-controls").hidden'), false);

  await b.evaluate('document.querySelector(\'[data-preset="gaussian"]\').click()');
  await b.evaluate("window.beforeShift = lesson.result.full.power.slice(); window.unshifted = lesson.result.full.image.slice()");
  await b.change("shift_x_a", 13);
  assert.ok(await b.evaluate("lesson.result.full.power.every((p,i) => Math.abs(p - beforeShift[i]) < 1e-14)"));
  assert.ok(await b.evaluate("lesson.result.full.image.some((p,i) => Math.abs(p - unshifted[i]) > 0.1)"));
  await b.change("remove_mean", true);
  assert.ok(await b.evaluate("Math.abs(lesson.result.full.real[0]) < 1e-14"));
  await b.evaluate('document.querySelector(\'[data-preset="noise"]\').click()');
  const peak = await b.evaluate("lesson.result.full.peak");
  await b.click("strongest");
  assert.equal(await b.value("mode_x"), String(peak.x)); assert.equal(await b.value("mode_y"), String(peak.y));
  await b.click("reset");
  await b.call("Emulation.setDeviceMetricsOverride", { width: 1440, height: 1800, deviceScaleFactor: 1, mobile: false });
  await b.navigate();
}).catch(error => { console.error(error); process.exitCode = 1; });
