"use strict";
const assert = require("node:assert/strict");
require("../TeachingTools/lesson_browser_test.js").runLessonBrowser("ModelEvidence", async b => {
  assert.ok(await b.hasColour("spectrum-plot", [213, 101, 52]));
  assert.ok(await b.hasColour("spectrum-plot", [80, 89, 166]));
  assert.ok(await b.hasColour("predictive-plot", [8, 127, 140]));
  assert.ok(await b.hasColour("predictive-plot", [80, 89, 166]));
  const flux = await b.evaluate("JSON.stringify(lesson.result.flux)");
  await b.change("prior", 0.1);
  const narrow = await b.text("log-bf"), estimate = await b.text("estimate"), maximum = await b.text("log-ml");
  await b.change("prior", 10);
  assert.equal(await b.evaluate("JSON.stringify(lesson.result.flux)"), flux);
  assert.equal(await b.text("estimate"), estimate);
  assert.equal(await b.text("log-ml"), maximum);
  assert.notEqual(await b.text("log-bf"), narrow);
  for (const values of [{ amplitude: -2, noise: 0.3 }, { amplitude: 3, noise: 2 }]) {
    await b.evaluate(`lesson.set(${JSON.stringify(values)})`); await b.healthy();
    assert.ok(await b.evaluate(`(() => { const r = lesson.result;
      return Number.isFinite(r.evidence.logBf) && r.evidence.probability >= 0 && r.evidence.probability <= 1 &&
        Math.abs(r.logEvidenceLine-r.logEvidenceNull-r.evidence.logBf)<1e-10;
    })()`));
    assert.ok(await b.hasColour("predictive-plot", [80, 89, 166]), "A broad predictive must remain visible.");
  }
}).catch(error => { console.error(error); process.exitCode = 1; });
