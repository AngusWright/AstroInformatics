"use strict";
const assert = require("node:assert/strict");
const { setTimeout: sleep } = require("node:timers/promises");
const { runBrowser } = require("../TeachingTools/browser_test.js");

async function suite({ evaluate, click, change, text, waitFor, healthy, navigate, openFile, screenshot, hasColour, call, base }) {
    const update = () => evaluate('document.getElementById("progress").value');
    assert.equal(await update(), 0);
    assert.equal(await text("playback-state"), "Paused");
    assert.ok(await hasColour("joint-plot", [199, 88, 0]));
    await click("step");
    assert.equal(await update(), 1);
    const firstDraw = await text("step-table");
    await click("sweep");
    assert.equal(await update(), 2);
    await click("rewind"); await click("step");
    assert.equal(await text("step-table"), firstDraw);
    await change("trail", 100); await change("burn", 20); await change("delay", 100);
    await click("rewind"); await click("step");
    assert.equal(await text("step-table"), firstDraw, "Display changes must not resample");
    await click("play");
    await waitFor('document.getElementById("progress").value >= 3');
    await click("play");
    const paused = await update();
    await sleep(180); assert.equal(await update(), paused);
    await change("sweeps", 10);
    assert.equal(await update(), 0);
    assert.equal(await evaluate('document.getElementById("burn").max'), "9");
    await change("jump_to", 19); await click("jump");
    await click("play"); assert.equal(await update(), 20);
    assert.equal(await text("play"), "Play");
    await click("play"); assert.equal(await update(), 1);
    await click("play");
    await change("jump_to", -1); await click("jump");
    assert.equal(await update(), 1);
    assert.equal(await evaluate('document.getElementById("error").hidden'), false);
    await click("rewind"); await healthy();

    await change("seed", -1);
    assert.equal(await evaluate('document.getElementById("play").disabled'), true);
    assert.equal(await evaluate('document.getElementById("error").hidden'), false);
    await change("seed", 666); await healthy();
    await change("sweeps", 100); await change("burn", 10);
    for (const family of ["gaussian", "student", "mixture", "poisson"]) {
      await change("family", family); await healthy();
      assert.equal(await update(), 0);
      assert.equal(await evaluate('document.getElementById("start_beta").matches(":disabled")'), family !== "poisson");
      await click("finish"); await click("tab-diagnostics"); await healthy();
      assert.match(await text("sample-count"), /100 completed sweeps; burn-in: 10; retained: 90/);
      for (const [j, colour] of [[0, [199, 88, 0]], [1, [0, 114, 178]]]) {
        for (const kind of ["trace", "acf"]) assert.ok(await hasColour(`${kind}-${j}`, colour), `${family}: ${kind} not drawn`);
        if (family !== "poisson") assert.ok(await hasColour(`marginal-${j}`, colour), `${family}: marginal overlay not drawn`);
      }
      await click("tab-explore"); await healthy();
      assert.ok(await hasColour("conditional-plot", [0, 114, 178]), `${family}: conditional not drawn`);
      if (family === "mixture") await screenshot("gibbs-desktop.png");
    }
    await change("start_beta", "");
    assert.equal(await evaluate('document.getElementById("play").disabled'), true);
    await change("family", "gaussian"); await healthy();
    for (const parameters of [
      { family: "mixture", rho: -0.98, separation: 8, weight: 0.05, scale1: 0.25, scale2: 3 },
      { family: "student", df: 3, rho: 0.98 },
      { family: "poisson", n: 100, total: 0, a: 1, hyper_shape: 1, hyper_rate: 0.1, start_beta: 1 }
    ]) {
      // Restore the deliberately invalid hidden start before selecting Poisson.
      if (parameters.family === "poisson") await change("start_beta", 1);
      for (const [id, value] of Object.entries(parameters)) await change(id, value);
      await click("finish"); await click("tab-diagnostics"); await healthy();
      await click("tab-explore"); await healthy();
    }
    await evaluate(`(() => {
      const tab = document.getElementById("tab-explore"); tab.focus();
      tab.dispatchEvent(new KeyboardEvent("keydown",{key:"ArrowRight",bubbles:true}));
    })()`);
    assert.equal(await evaluate("document.activeElement.id"), "tab-diagnostics");
    await click("tab-explore");
    await call("Emulation.setDeviceMetricsOverride", { width: 390, height: 844, deviceScaleFactor: 2, mobile: false });
    await navigate(base);
    assert.equal(await evaluate('document.getElementById("settings").open'), false);
    assert.ok(await evaluate("document.documentElement.scrollWidth <= innerWidth"), "Mobile viewport overflows");
    await click("finish"); await screenshot("gibbs-mobile.png"); await healthy();

    await openFile();
    await click("step"); assert.equal(await update(), 1); await healthy();
    console.log("Browser passed: four families, Canvas plots, timers, replay, burn-in, invalid inputs, keyboard tabs, mobile layout, offline file use, no runtime downloads.");
}
runBrowser("Gibbs", suite).catch(error => { console.error(error); process.exitCode = 1; });
