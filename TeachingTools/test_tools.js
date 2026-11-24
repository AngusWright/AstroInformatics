"use strict";
const fs = require("node:fs");
const path = require("node:path");
const { spawnSync } = require("node:child_process");

const ROOT = path.resolve(__dirname, "..");
const tools = JSON.parse(fs.readFileSync(path.join(ROOT, "html-tools.json"), "utf8"));
const args = process.argv.slice(2);
if (args.length > 1 || (args.length === 1 && args[0] !== "--browser")) {
  throw new Error("Usage: node TeachingTools/test_tools.js [--browser]");
}
const browser = args[0] === "--browser";
const tests = tools.flatMap(tool => {
  const directory = path.join(ROOT, tool.path);
  const selected = browser ? ["test_browser.js"] :
    fs.readdirSync(directory).filter(name => /^test_.*\.js$/.test(name) && name !== "test_browser.js").sort();
  if (!selected.length || selected.some(name => !fs.statSync(path.join(directory, name)).isFile())) {
    throw new Error(`${tool.path}: missing ${browser ? "browser" : "numerical"} regression tests`);
  }
  return selected.map(name => `${tool.path}/${name}`);
});

function run(arguments_) {
  const child = spawnSync(process.execPath, arguments_, { cwd: ROOT, stdio: "inherit" });
  if (child.error) throw child.error;
  if (!Number.isInteger(child.status)) throw new Error(`Test process terminated by ${child.signal}`);
  if (child.status !== 0) process.exit(child.status);
}
if (browser) {
  for (const file of tests) run([file]);
  run(["TeachingTools/test_gallery_browser.js"]);
} else {
  run(["--test", "--test-reporter=spec", "TeachingTools/test_statistics.js", "TeachingTools/test_export.js", ...tests]);
}
