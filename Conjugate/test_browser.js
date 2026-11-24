"use strict";
const assert = require("node:assert/strict");
const { runBrowser } = require("../TeachingTools/browser_test.js");

runBrowser("Conjugate", async b => {
  const { click, change, text, evaluate, healthy } = b;
  assert.equal(await text("sample-size"), "10");
  assert.ok(await b.hasColour("density-plot", [199, 88, 0]));
  const data = await text("data-description");
  await change("prior_mean", 4);
  assert.equal(await text("data-description"), data, "A prior change must not change observations");
  await click("strong_prior");
  assert.ok(parseFloat(await text("prior-weight")) > 98);
  await click("many_data");
  assert.equal(await text("sample-size"), "10,000");
  assert.ok(parseFloat(await text("prior-weight")) < 1);
  await click("weak_prior");
  assert.ok(parseFloat(await text("prior-weight")) < 1);
  await click("reset");
  for (const model of ["normal", "binomial", "poisson", "gamma"]) {
    await change("model", model); await healthy();
    await change("no_data", true);
    assert.equal(await text("sample-size"), "0");
    const rows = await evaluate('[...document.querySelectorAll("#summary tr")].map(r=>[...r.cells].slice(1).map(c=>c.textContent))');
    assert.deepEqual(rows[0], rows[1]);
    assert.equal(await evaluate('document.getElementById("likelihood-legend").hidden'), true);
    await change("no_data", false);
    await click("tab-role"); await healthy();
    assert.ok(await b.hasColour("role-plot", [199, 88, 0]));
    assert.ok(await b.hasColour("weight-plot", [0, 114, 178]));
    await click("tab-update");
  }
  for (const parameters of [
    { model: "binomial", log_n: 4, beta_a: 0.5, beta_b: 0.5, success_fraction: 0 },
    { model: "binomial", success_fraction: 1 },
    { model: "poisson", gamma_a: 0.5, gamma_b: 0.1, count_mean: 0 },
    { model: "gamma", shape: 10, gamma_mean: 0.1, gamma_b: 20 },
    { model: "normal", sigma: 0.1, prior_sd: 0.05, prior_mean: -5, sample_mean: 5 }
  ]) {
    for (const [id, value] of Object.entries(parameters)) await change(id, value);
    await healthy(); await click("tab-role"); await healthy(); await click("tab-update");
  }
  await change("model", "unknown");
  assert.equal(await evaluate('document.getElementById("error").hidden'), false);
  assert.equal(await evaluate('document.querySelectorAll("#summary tr").length'), 0);
  await change("model", "normal"); await healthy();
  await click("reset"); await b.screenshot("conjugate-desktop.png");
  await click("tab-explain");
  const sufficiency = await evaluate('[...document.querySelectorAll("#sufficiency tr")].map(r=>[...r.cells].slice(1).map(c=>c.textContent))');
  assert.deepEqual(sufficiency[0], sufficiency[1]);
  await click("tab-update");
  await evaluate(`(() => {
    const tab=document.getElementById("tab-update"); tab.focus();
    tab.dispatchEvent(new KeyboardEvent("keydown",{key:"ArrowRight",bubbles:true}));
  })()`);
  assert.equal(await evaluate("document.activeElement.id"), "tab-role");
  await b.call("Emulation.setDeviceMetricsOverride", { width: 390, height: 844, deviceScaleFactor: 2, mobile: false });
  await b.navigate();
  assert.equal(await evaluate('document.getElementById("settings").open'), false);
  assert.ok(await evaluate("document.documentElement.scrollWidth <= innerWidth"));
  await b.screenshot("conjugate-mobile.png");
  await b.openFile(); await change("model", "gamma"); await healthy();
  await b.call("Page.navigate", { url: "about:blank" });
  await b.waitFor('document.URL === "about:blank" && document.readyState === "complete"');
  await evaluate("document.documentElement.remove()");
  await b.navigate();
  assert.equal(await text("sample-size"), "10");
}).catch(error => { console.error(error); process.exitCode = 1; });
