"use strict";
const assert = require("node:assert/strict");
const { runLessonBrowser } = require("../TeachingTools/lesson_browser_test.js");

runLessonBrowser("InverseSampling", async b => {
  assert.equal(await b.evaluate("document.title"), "Turning Uniform Randomness into a Universe");
  assert.equal(await b.evaluate('document.querySelector(".lecture-link").href'),
    "https://anguswright.github.io/LectureNotes/IntroductionToStatistics_Section2d_html.html");
  assert.equal(await b.text("revealed"), "500");
  assert.ok(await b.hasColour("cdf-plot", [80, 89, 166]));
  assert.ok(await b.hasColour("density-plot", [8, 127, 140]));
  await b.evaluate("window.originalUniforms = lesson.result.data");
  const initialImage = await b.evaluate('document.getElementById("density-plot").toDataURL()');
  await b.change("slope", 1);
  assert.equal(await b.evaluate("lesson.result.data === originalUniforms"), true);
  assert.equal(await b.text("target-mass"), "1.000000");
  assert.notEqual(await b.evaluate('document.getElementById("density-plot").toDataURL()'), initialImage);
  assert.ok(await b.evaluate(`Math.abs(InverseSamplingModel.logDensity(lesson.result.target, -0.5) -
    InverseSamplingModel.logDensity(lesson.result.target, 1)) < 1e-12`));
  for (const slope of [0.5, 3.5]) { await b.change("slope", slope); await b.healthy(); }
  for (const family of ["exponential", "lognormal"]) {
    await b.change("family", family); await b.healthy();
    assert.equal(await b.evaluate("lesson.result.data === originalUniforms"), true);
    assert.ok(Number(await b.text("target-mass")) < 1);
    assert.ok(Number(await b.text("target-mass")) > 0.9999);
    assert.match(await b.text("target-description"), family === "exponential" ? /Exponential/ : /Log-normal/);
    if (family === "exponential") {
      const [u, x] = await b.evaluate("[lesson.result.uniforms[0], lesson.result.values[0]]");
      assert.ok(Math.abs(x + Math.log1p(-u)) < 1e-12);
      const values = await b.evaluate("lesson.result.values");
      await b.change("slope", 0.5);
      assert.deepEqual(await b.evaluate("lesson.result.values"), values, "Inactive power exponent must not change exponential draws");
    }
  }
  for (const n of [50, 1500]) {
    await b.change("n", n); await b.healthy();
    assert.equal(await b.text("revealed"), String(n));
    assert.equal(await b.evaluate("lesson.result.data === originalUniforms"), true);
  }
  await b.change("n", 100);
  await b.click("rewind");
  assert.equal(await b.text("revealed"), "0");
  assert.equal(await b.text("ks"), "No draws yet");
  assert.ok(await b.hasColour("density-plot", [8, 127, 140]), "The analytic target remains visible without observations");
  assert.equal(await b.hasColour("density-plot", [233, 171, 142]), false, "No histogram bars before a draw");
  assert.equal(await b.hasColour("cdf-plot", [8, 127, 140]), false, "No transformed dots before a draw");
  await b.click("step");
  assert.equal(await b.text("revealed"), "1");
  assert.ok(await b.hasColour("cdf-plot", [8, 127, 140]));
  const firstValue = await b.evaluate("lesson.result.values[0]");
  await b.evaluate('document.getElementById("frame").value = "50"; document.getElementById("frame").dispatchEvent(new Event("input", { bubbles:true }))');
  assert.equal(await b.text("revealed"), "50");
  await b.click("finish");
  assert.equal(await b.evaluate("lesson.result.values[0]"), firstValue);
  assert.equal(await b.evaluate("lesson.result.data === originalUniforms"), true);
  await b.click("redraw");
  await b.evaluate("window.newUniforms = lesson.result.data");
  await b.change("family", "power");
  assert.equal(await b.evaluate("lesson.result.data === newUniforms && lesson.parameters.seed === 43"), true);
  await b.change("seed", 0); await b.healthy();
  await b.change("seed", 2147483647); await b.healthy();
  await b.change("family", "unknown");
  assert.equal(await b.evaluate('document.getElementById("error").hidden'), false);
  await b.click("reset"); await b.healthy();
}).catch(error => { console.error(error); process.exitCode = 1; });
