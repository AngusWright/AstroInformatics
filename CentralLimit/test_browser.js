"use strict";
const assert = require("node:assert/strict");
const { runLessonBrowser } = require("../TeachingTools/lesson_browser_test.js");

runLessonBrowser("CentralLimit", async b => {
  assert.equal(await b.evaluate("document.title"), "Why Stacking Works (and When It Doesn't)");
  assert.equal(await b.evaluate('document.querySelector(".lecture-link").href'),
    "https://anguswright.github.io/LectureNotes/IntroductionToStatistics_Section2c_html.html");
  assert.equal(await b.text("stack-n"), "30");
  assert.ok(await b.hasColour("population-plot", [8, 127, 140]));
  assert.ok(await b.hasColour("sampling-plot", [80, 89, 166]));
  const initialIqr = Number(await b.text("empirical-iqr"));
  const populationImage = await b.evaluate('document.getElementById("population-plot").toDataURL()');
  const samplingImage = await b.evaluate('document.getElementById("sampling-plot").toDataURL()');
  await b.evaluate("window.originalStackData = lesson.result.data");
  await b.change("n", 100);
  assert.equal(await b.evaluate("lesson.result.data === originalStackData"), true);
  assert.equal(await b.evaluate('document.getElementById("population-plot").toDataURL()'), populationImage);
  assert.notEqual(await b.evaluate('document.getElementById("sampling-plot").toDataURL()'), samplingImage);
  assert.ok(Number(await b.text("empirical-iqr")) < initialIqr);
  assert.equal(await b.text("mean-sd"), "0.1");
  await b.change("repetitions", 200);
  assert.equal(await b.text("experiments"), "200");
  assert.equal(await b.evaluate("lesson.result.data === originalStackData"), true);
  await b.change("repetitions", 3000);
  for (const population of ["mixture", "exponential", "cauchy"]) {
    await b.change("population", population); await b.healthy();
    assert.ok(await b.hasColour("sampling-plot", [80, 89, 166]));
    assert.equal(await b.evaluate("CentralLimitModel.atFrame(lesson.result).finite"), population !== "cauchy");
  }
  assert.equal(await b.text("mean-sd"), "Undefined");
  assert.match(await b.text("moment-note"), /mean and variance are undefined/);
  assert.match(await b.text("reference-label"), /exact standard Cauchy/);
  const cauchyWidth = Number(await b.text("empirical-iqr"));
  await b.evaluate("window.cauchyStackData = lesson.result.data");
  await b.click("rewind");
  assert.equal(await b.text("stack-n"), "1");
  assert.ok(cauchyWidth / Number(await b.text("empirical-iqr")) > 0.7);
  await b.click("finish");
  assert.equal(Number(await b.text("empirical-iqr")), cauchyWidth);
  assert.equal(await b.evaluate("lesson.result.data === cauchyStackData"), true);
  await b.click("redraw");
  await b.evaluate("window.newStackData = lesson.result.data");
  await b.change("n", 1);
  assert.equal(await b.evaluate("lesson.result.data === newStackData && lesson.parameters.seed === 43"), true);
  await b.change("seed", 0); await b.healthy();
  await b.change("seed", 2147483647); await b.healthy();
  await b.change("population", "unknown");
  assert.equal(await b.evaluate('document.getElementById("error").hidden'), false);
  await b.click("reset"); await b.healthy();
}).catch(error => { console.error(error); process.exitCode = 1; });
