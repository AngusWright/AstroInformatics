"use strict";
const assert = require("node:assert/strict");
const { runLessonBrowser } = require("../TeachingTools/lesson_browser_test.js");

runLessonBrowser("DensityEstimation", async b => {
  assert.equal(await b.evaluate("document.title"), "Can You Trust That Bump?");
  assert.equal(await b.evaluate('document.querySelector(".lecture-link").href'),
    "https://anguswright.github.io/LectureNotes/IntroductionToStatistics_Section1b_html.html");
  assert.equal(await b.text("objects"), "200"); assert.equal(await b.text("hist-area"), "1.000");
  assert.ok(await b.hasColour("density-plot", [8, 127, 140]));
  assert.ok(await b.hasColour("density-plot", [80, 89, 166]));
  assert.ok(await b.hasColour("density-plot", [203, 216, 235]));
  await b.evaluate("window.originalColourData = lesson.result.data");
  const values = await b.evaluate("lesson.result.values");
  const density = await b.evaluate("lesson.result.density");
  const counts = await b.evaluate("lesson.result.histogram.counts");
  const image = await b.evaluate('document.getElementById("density-plot").toDataURL()');
  await b.change("origin", 0.75);
  assert.equal(await b.evaluate("lesson.result.data === originalColourData"), true);
  assert.deepEqual(await b.evaluate("lesson.result.values"), values);
  assert.deepEqual(await b.evaluate("lesson.result.density"), density);
  assert.notDeepEqual(await b.evaluate("lesson.result.histogram.counts"), counts);
  assert.notEqual(await b.evaluate('document.getElementById("density-plot").toDataURL()'), image);
  await b.change("origin", 0);
  const originZero = await b.evaluate("lesson.result.histogram.counts");
  await b.change("origin", 1);
  assert.deepEqual(await b.evaluate("lesson.result.histogram.counts"), originZero);
  for (const [id, value] of [["width", 0.1], ["width", 1.5], ["bandwidth", 0.05], ["bandwidth", 1]]) {
    await b.change(id, value); await b.healthy();
    assert.deepEqual(await b.evaluate("lesson.result.values"), values);
    assert.equal(await b.evaluate("lesson.result.data === originalColourData"), true);
    assert.equal(await b.text("hist-area"), "1.000");
  }
  assert.ok(Number(await b.text("outside-mass")) > 1, "The broad KDE has real probability beyond the viewport");
  await b.evaluate('document.querySelector(\'[data-preset="under"]\').click()');
  const sharpPeaks = Number(await b.text("peaks"));
  await b.evaluate('document.querySelector(\'[data-preset="over"]\').click()');
  assert.equal(await b.text("peaks"), "1");
  assert.ok(sharpPeaks > 1);
  assert.deepEqual(await b.evaluate("lesson.result.values"), values);
  await b.change("n", 50); assert.equal(await b.text("objects"), "50");
  assert.deepEqual(await b.evaluate("lesson.result.values"), values.slice(0, 50));
  await b.change("n", 1000); assert.equal(await b.text("objects"), "1000");
  assert.deepEqual(await b.evaluate("lesson.result.values.slice(0, 200)"), values);
  await b.click("redraw");
  await b.evaluate("window.newColourData = lesson.result.data");
  await b.change("bandwidth", 0.05);
  assert.equal(await b.evaluate("lesson.result.data === newColourData && lesson.parameters.seed === 43"), true);
  await b.change("seed", 0); await b.healthy();
  await b.change("seed", 2147483647); await b.healthy();
  await b.change("seed", 1.5);
  assert.equal(await b.evaluate('document.getElementById("error").hidden'), false);
  assert.equal(await b.text("objects"), "—");
  await b.click("reset"); await b.healthy();
}).catch(error => { console.error(error); process.exitCode = 1; });
