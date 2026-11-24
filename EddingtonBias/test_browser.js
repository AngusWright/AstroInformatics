"use strict";
const assert = require("node:assert/strict");
require("../TeachingTools/lesson_browser_test.js").runLessonBrowser("EddingtonBias", async b => {
  assert.ok(await b.hasColour("mass-plot", [8, 127, 140]));
  assert.ok(await b.hasColour("mass-plot", [213, 101, 52]));
  await b.change("noise", 0);
  assert.equal(await b.text("up"), "0"); assert.equal(await b.text("down"), "0");
  assert.ok(await b.evaluate("lesson.result.truth.every((x,i)=>x===lesson.result.final.observed[i])"));
  const truth = await b.evaluate("JSON.stringify(lesson.result.truth)");
  await b.change("noise", 0.6);
  assert.equal(await b.evaluate("JSON.stringify(lesson.result.truth)"), truth);
  const measured = await b.evaluate("JSON.stringify(lesson.result.final.observed)");
  await b.change("threshold", 11.5);
  assert.equal(await b.evaluate("JSON.stringify(lesson.result.final.observed)"), measured);
  assert.ok(await b.evaluate(`(() => { const s = lesson.result.final;
    return s.observedHigh === s.trueHigh + s.up - s.down &&
      s.lowLow + s.up + s.down + s.highHigh === 8000 &&
      s.histogram.excluded === 0;
  })()`));
  await b.click("rewind");
  assert.equal(await b.text("up"), "0"); assert.equal(await b.text("down"), "0");
  await b.click("finish");
  assert.equal(await b.evaluate("JSON.stringify(lesson.result.final.observed)"), measured);
  for (const values of [{ slope: 0.3, threshold: 9.5 }, { slope: 3, threshold: 11.5 }]) {
    await b.evaluate(`lesson.set(${JSON.stringify(values)})`); await b.healthy();
    assert.equal(await b.evaluate("lesson.result.final.histogram.counts.reduce((a,b)=>a+b,0)"), 8000);
    assert.ok(await b.hasColour("mass-plot", [213, 101, 52]));
  }
}).catch(error => { console.error(error); process.exitCode = 1; });
