"use strict";
const assert = require("node:assert/strict");
const { runBrowser } = require("../TeachingTools/browser_test.js");

runBrowser("ShiftFitting", async b => {
  await b.healthy();
  assert.equal(await b.evaluate("document.title"), "Shift Fitting under Heteroskedasticity");
  assert.match(await b.text("fit-table"), /±/);
  assert.equal(await b.evaluate("document.querySelectorAll('#fit-table tr').length"), 3);
  const observations = await b.evaluate("document.getElementById('main').dataset.observations");
  const generation = await b.evaluate("document.getElementById('main').dataset.generation");
  const initialTable = await b.text("fit-table");

  for (const mode of ["bias", "correlation", "data_model"]) {
    await b.change("plot_type", mode);
    await b.healthy();
    assert.equal(await b.evaluate("document.getElementById('main').dataset.observations"), observations);
  }
  await b.change("realisation", 2);
  assert.notEqual(await b.text("fit-table"), initialTable);
  await b.change("annot_cex", 2);
  for (const weighting of ["equal", "known", "legacy"]) {
    await b.change("weighting", weighting);
    await b.healthy();
    assert.equal(await b.evaluate("document.getElementById('main').dataset.observations"), observations);
  }
  await b.change("model_min", .8);
  await b.change("nmodels", 4);
  await b.healthy();
  assert.equal(await b.evaluate("document.querySelectorAll('#fit-table tr').length"), 4);
  assert.equal(await b.evaluate("document.getElementById('main').dataset.generation"), generation);
  await b.click("generate");
  assert.equal(await b.evaluate("document.getElementById('main').dataset.observations"), observations, "Same seed must exactly reproduce observations");

  // Data controls are staged, never silently applied to an existing experiment.
  await b.change("seed", 2026);
  await b.change("model_type", "linear");
  assert.equal(await b.evaluate("document.getElementById('main').dataset.observations"), observations);
  assert.match(await b.text("pending"), /Click Generate/);
  await b.change("n_noise", 1);
  await b.change("num_samples", 1.5);
  await b.change("model_min", .9);
  await b.change("nmodels", 3);
  for (const noise of ["constant", "linear", "log-linear"]) {
    await b.change("noise_model", noise);
    await b.click("generate");
    await b.healthy();
    assert.match(await b.text("truth-summary"), /intercept/);
    assert.match(await b.text("fit-table"), /Not fitted/);
    for (const mode of ["bias", "correlation", "data_model"]) {
      await b.change("plot_type", mode);
      await b.healthy();
    }
  }
  await b.change("noise_model", "constant");
  await b.click("generate");
  await b.change("correct_htrsk", true);
  await b.healthy();
  assert.match(await b.text("interval-note"), /no analytic|noisy fits/i);
  assert.match(await b.text("fit-table"), /Unavailable \(n=0\)/);
  await b.change("correct_htrsk", false);

  // No stale fit is presented as the answer to an invalid input.
  await b.change("model_min", "");
  assert.equal(await b.evaluate("document.getElementById('error').hidden"), false);
  assert.equal(await b.text("fit-table"), "");
  await b.change("plot_type", "bias");
  assert.equal(await b.evaluate("document.getElementById('error').hidden"), false, "A display change must not hide an unresolved inference error");
  await b.change("model_min", .9);
  await b.healthy();
  await b.change("sample_min", 5);
  await b.click("generate");
  assert.match(await b.text("error"), /sampling range/);
  await b.change("sample_min", -4);
  await b.click("generate");
  await b.healthy();
  await b.change("seed", -1);
  await b.click("generate");
  assert.match(await b.text("error"), /Seed/);
  await b.change("seed", 2026);
  await b.change("n_noise", 0);
  await b.click("generate");
  await b.healthy();
  await b.change("plot_type", "correlation");
  assert.match(await b.text("plot-note"), /need at least two/);
  await b.change("model_type", "gaussian");
  await b.change("num_samples", 0);
  await b.click("generate");
  assert.match(await b.text("error"), /cannot identify/);
  await b.change("num_samples", 1.5);
  await b.change("n_noise", 1);
  await b.click("generate");
  await b.healthy();
  await b.change("noise_model", "constant");
  await b.change("seed", 666);
  await b.click("generate");
  await b.change("correct_htrsk", true);
  await b.healthy();
  assert.match(await b.text("fit-table"), /Unavailable \(n=0\)/);
  await b.change("correct_htrsk", false);
  await b.change("plot_type", "data_model");
  assert.equal(await b.hasColour("fit-plot", [24, 126, 145]), true);
  await b.screenshot("shift-fitting-desktop.png");

  await b.call("Emulation.setDeviceMetricsOverride", { width: 390, height: 844, deviceScaleFactor: 2, mobile: true });
  await b.waitFor("document.getElementById('fit-plot').width >= document.getElementById('fit-plot').clientWidth * 1.9");
  assert.equal(await b.evaluate("document.documentElement.scrollWidth <= window.innerWidth + 1"), true);
  assert.equal(await b.evaluate("document.getElementById('fit-plot').clientWidth <= window.innerWidth"), true);
  await b.healthy();
  await b.navigate();
  assert.equal(await b.evaluate('document.getElementById("settings").open'), false);
  await b.screenshot("shift-fitting-mobile.png");
  await b.openFile();
  await b.healthy();
  assert.match(await b.text("fit-table"), /±/);
  await b.change("weighting", "known");
  await b.healthy();
});
