"use strict";
const assert = require("node:assert/strict");
require("../TeachingTools/lesson_browser_test.js").runLessonBrowser("Correlation", async b => {
  assert.equal(await b.evaluate("document.title"), "Confounding Correlations");
  const fullR = await b.text("full-r");
  assert.ok(Number(fullR) > 0.6);
  assert.ok(Math.abs(Number(await b.text("intrinsic-r"))) < 0.08);
  assert.ok(await b.hasColour("flux-plot", [80, 89, 166]), "The selected-slice regression must be drawn");
  assert.ok(await b.evaluate(`(() => {
    const c = document.getElementById("intrinsic-plot"), p = c.getContext("2d").getImageData(0,0,c.width,c.height).data;
    let coloured = 0;
    for (let i=0;i<p.length;i+=4) if (Math.max(p[i],p[i+1],p[i+2])-Math.min(p[i],p[i+1],p[i+2])>30) coloured++;
    return coloured > 200;
  })()`), "The intrinsic plot must contain the actual galaxy points");
  await b.evaluate(`globalThis.populationBase = lesson.result.base;
    globalThis.observedX = JSON.stringify(lesson.result.x);
    globalThis.observedY = JSON.stringify(lesson.result.y);
    Math.random = () => { throw new Error("Unseeded randomness"); };`);
  await b.change("width", 0.05); await b.healthy();
  assert.equal(await b.text("full-r"), fullR);
  assert.equal(await b.evaluate("lesson.result.base === populationBase && JSON.stringify(lesson.result.x) === observedX && JSON.stringify(lesson.result.y) === observedY"), true);
  assert.ok(await b.evaluate("lesson.result.slice.n < 150 && Math.abs(lesson.result.slice.correlation) < lesson.result.full.correlation"));
  await b.change("centre", 0); await b.change("width", 1);
  assert.equal(await b.evaluate("lesson.result.lower === 0 && lesson.result.upper === 0.5 && lesson.result.slice.n < 1800"), true);
  await b.change("scatter", 0); await b.healthy();
  assert.equal(await b.text("full-r"), "1.00");
  assert.equal(await b.text("intrinsic-r"), "Undefined");
  assert.match(await b.text("variance-note"), /zero variance/);
  await b.change("scatter", 2); await b.change("centre", 1); await b.change("width", 0.05);
  await b.change("seed", 2147483647); await b.healthy();
  assert.ok(await b.evaluate("lesson.result.slice.n >= 3 && Number.isFinite(lesson.result.slice.correlation)"));
  await b.click("reset"); await b.click("redraw");
  await b.evaluate("globalThis.newPopulation = lesson.result.base");
  await b.change("width", 0.1);
  assert.equal(await b.evaluate("lesson.result.base === newPopulation && lesson.generation === 1"), true,
    "Slicing after a new realisation must not silently return to the base seed");
  await b.click("reset");
}).catch(error => { console.error(error); process.exitCode = 1; });
