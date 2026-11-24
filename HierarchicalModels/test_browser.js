"use strict";
const assert = require("node:assert/strict");
require("../TeachingTools/lesson_browser_test.js").runLessonBrowser("HierarchicalModels", async b => {
  assert.ok(await b.hasColour("stars-plot", [80, 89, 166]));
  assert.ok(await b.hasColour("stars-plot", [8, 127, 140]));
  const observed = await b.evaluate("JSON.stringify(lesson.result.observed)");
  await b.change("pooling", "none");
  assert.equal(await b.text("shift"), "0.00");
  assert.equal(await b.text("population"), "Not inferred");
  assert.equal(await b.text("weight"), "100.0%");
  assert.equal(await b.evaluate("JSON.stringify(lesson.result.observed)"), observed);
  await b.change("pooling", "complete");
  assert.equal(await b.evaluate("JSON.stringify(lesson.result.observed)"), observed);
  assert.ok(await b.evaluate(`lesson.result.posterior.mean.every(x=>x===lesson.result.posterior.populationMean) &&
    lesson.result.posterior.variance.every(x=>x>0)`));
  await b.change("scatter", 0);
  const complete = await b.evaluate("JSON.stringify(lesson.result.posterior)");
  await b.change("pooling", "partial");
  assert.equal(await b.evaluate("JSON.stringify(lesson.result.posterior)"), complete);
  assert.equal(await b.text("weight"), "0.0%");
  await b.click("rewind");
  assert.equal(await b.text("shift"), "0.00");
  assert.match(await b.text("phase-note"), /not posterior credible intervals/);
  await b.click("finish");
  for (const noise of [0.3, 3]) {
    await b.evaluate(`lesson.set({scatter:3,noise:${noise}})`); await b.healthy();
    assert.ok(await b.evaluate("lesson.result.posterior.variance.every(x=>Number.isFinite(x)&&x>0)"));
    assert.ok(await b.hasColour("stars-plot", [80, 89, 166]));
  }
}).catch(error => { console.error(error); process.exitCode = 1; });
