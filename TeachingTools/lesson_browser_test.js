"use strict";
const assert = require("node:assert/strict");
const { runBrowser } = require("./browser_test.js");

function runLessonBrowser(page, verify) {
  return runBrowser(page, async b => {
    const metrics = () => b.evaluate('[...document.querySelectorAll("[data-result]")].map(node => node.textContent)');
    const original = await metrics();
    assert.ok(original.length > 0 && original.every(value => value.trim().length > 0), "A lesson must show numerical results");
    assert.equal(await b.evaluate('document.getElementById("play")?.getAttribute("aria-pressed") || "false"'), "false");
    if (verify) await verify(b);
    for (const preset of await b.evaluate('[...document.querySelectorAll("[data-preset]")].map(node => node.dataset.preset)')) {
      await b.evaluate(`[...document.querySelectorAll("[data-preset]")].find(node =>
        node.dataset.preset === ${JSON.stringify(preset)}).click()`);
      await b.healthy();
    }
    await b.click("reset");
    assert.deepEqual(await metrics(), original, "Reset must restore seeded results");
    await b.click("redraw"); await b.healthy(); await b.click("reset");
    assert.deepEqual(await metrics(), original, "New draws must not consume global randomness or change the reset result");
    if (await b.evaluate('!!document.getElementById("play")')) {
      await b.click("rewind"); assert.equal(await b.value("frame"), "0");
      await b.click("step"); assert.equal(await b.value("frame"), "1");
      await b.click("play"); assert.ok(Number(await b.value("frame")) >= 2);
      await b.click("play");
      const paused = await b.value("frame");
      await b.evaluate("new Promise(resolve => setTimeout(resolve, 250))");
      assert.equal(await b.value("frame"), paused, "Pause must stop the timer");
      await b.click("finish"); await b.healthy();
      assert.deepEqual(await metrics(), original, "Playback must reveal cached results, not resample them");
    }
    await b.change("seed", -1);
    assert.equal(await b.evaluate('document.getElementById("error").hidden'), false);
    assert.ok((await metrics()).every(value => value === "\u2014"), "Invalid inputs must not leave stale results");
    await b.click("reset"); await b.healthy();
    await b.click("tab-explain");
    assert.equal(await b.evaluate('document.getElementById("panel-explain").hidden'), false);
    await b.evaluate(`(() => {
      const tab = document.getElementById("tab-explain"); tab.focus();
      tab.dispatchEvent(new KeyboardEvent("keydown", { key: "ArrowRight", bubbles: true }));
    })()`);
    assert.equal(await b.evaluate("document.activeElement.id"), "tab-explore");
    await b.screenshot(`${page.toLowerCase()}-desktop.png`);
    await b.call("Emulation.setDeviceMetricsOverride", { width: 390, height: 844, deviceScaleFactor: 2, mobile: false });
    await b.navigate();
    assert.equal(await b.evaluate('document.getElementById("settings").open'), false);
    assert.ok(await b.evaluate("document.documentElement.scrollWidth <= innerWidth"), "Mobile layout must not overflow");
    await b.screenshot(`${page.toLowerCase()}-mobile.png`);
    await b.openFile(); await b.healthy();
  });
}
module.exports = { runLessonBrowser };
