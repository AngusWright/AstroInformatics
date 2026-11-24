"use strict";
const assert = require("node:assert/strict");
const { runBrowser } = require("./browser_test.js");
const tools = require("../html-tools.json");

process.env.HTML_SITE_ROOT = process.env.HTML_SITE_ROOT || "_site";
runBrowser("", async b => {
  const cards = await b.evaluate('[...document.querySelectorAll(".card h3 a")].map(link => [link.textContent, link.getAttribute("href")])');
  assert.deepEqual(cards, tools.map(tool => [tool.title, `${tool.path}/`]));
  const headings = await b.evaluate('[...document.querySelectorAll("main section > h2")].map(node => node.textContent)');
  assert.deepEqual(headings, [...new Set(tools.map(tool => tool.group))]);
  assert.ok(await b.evaluate("document.documentElement.scrollWidth <= innerWidth"));
  await b.screenshot("gallery-desktop.png");
  await b.call("Emulation.setDeviceMetricsOverride", { width: 390, height: 844, deviceScaleFactor: 2, mobile: false });
  assert.ok(await b.evaluate("document.documentElement.scrollWidth <= innerWidth"), "Gallery must not overflow on mobile");
  assert.ok(await b.evaluate('[...document.querySelectorAll(".card")].every(node => node.getBoundingClientRect().width <= innerWidth)'));
  await b.screenshot("gallery-mobile.png");
  await b.openFile();
  assert.equal(await b.evaluate('document.querySelectorAll(".card").length'), tools.length);
}).catch(error => { console.error(error); process.exitCode = 1; });
